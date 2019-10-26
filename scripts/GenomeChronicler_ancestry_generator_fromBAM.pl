#!/usr/bin/env perl -w

package GC_Ancestry;


use warnings;
use strict;
use File::Basename;


use Exporter qw(import);
 
#our @EXPORT_OK = qw(add multiply);



my $dir="/GenomeChronicler/";
my $resultsdir = $dir;

#- Also needs tabix and R in the PATH

my $plink="${dir}/software/plink";
my $gatk="${dir}/software/GenomeAnalysisTK.jar";
my $ref_hs38="${dir}/reference/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa";
my $initialBIM="${dir}/reference/1kGP_GRCh38_exome.bim";
my $initialAnc="${dir}/reference/1kGP_GRCh38_exome";

my $bcftools="bcftools";
my $samtools="samtools";
my $bgzip="bgzip";
my $tabix="tabix";

my $hasCHRflag=0;
my $numThreads = 4;


if(!@ARGV) {
	die "Please provide a single base recalibrated BAM file to this script\n";
}

my $bam=$ARGV[0];
(my $sample= basename($bam)) =~ s/\.[^.]+$//;
$sample =~ s/\.recal//g;
$sample =~ s/\.bam\.clean//gi;

$resultsdir = $ARGV[1] if(defined($ARGV[1]));

$numThreads = $ARGV[2] if(defined($ARGV[2]));

system("mkdir -p ${resultsdir}/results/results_${sample}/temp");

#Also check for the need to samtools index this stuff
#java -jar software/picard.jar CreateSequenceDictionary R=software/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa O=software/GRCh38_full_analysis_set_plus_decoy_hla_noChr.dict


#------------------------------------------
#-- I had already preprocessed the 1kGP data to reduce it to a PLINK file with
#--	only SNPs from exome (914577 total)
#-- But from these, here I remove all the A/T or C/G SNPs to make sure that there
#--	will not be artifacts due to reporting SNPS i different stransd ib both datasets
#-- I also remove sexual chromosomes
#------------------------------------------

&_get_list_of_snps_from_1kGP();

#------------------------------------------
#-- Convert BAM to a gvcf that contains only
#-  	the 1kgp snps filtered above
#------------------------------------------

&_get_gvcf();

#------------------------------------------
#-- Convert GVCF data to PLINK format
#------------------------------------------

&_gvcf_to_plink();

#------------------------------------------
#-- Merge gvcf with 1kGP data. I need to merge the files twice, the first one will fail if there are SNPs with 
#--	different alleles in both datasets. The second merge will exclude these SNPs 
#------------------------------------------

&_merge_pgp_1kGP();

#------------------------------------------
#-- Reduce the merged dataset further by removing SPs in LD
#------------------------------------------

&_subset_unlinked_snps();


#------------------------------------------
#-- Perform PCA
#------------------------------------------

my $runstr6="$plink";
$runstr6.=" --bfile ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP_pruned";
$runstr6.=" --out ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP_pruned_pca_20";
$runstr6.=" --pca";
`$runstr6`;



#############################################################################
###########                   SUBROUTINES                            ########
#############################################################################


sub _get_list_of_snps_from_1kGP {

	#---------------------------------------
	#- get a list of SNPs from 1kGP
	#-	(no A/t, C/G or sex SNPs)
	#---------------------------------------


	open(INF, "${initialBIM}");
	open(OUTF, ">${resultsdir}/results/results_${sample}/temp/1kGP_GRCh38_nonAT_CG.bed");
	print OUTF "#CHROM\tST0\tPOS\tID\n";
	while(<INF>){
		chomp $_;
	
		my @a=split("\t",$_);


		if($hasCHRflag == 1) {
			if($a[0] !~ m/^chr/) {
				$a[0] = "chr".$a[0];
			}
		}
		
		my $pp = $a[0];
		$pp =~ s/chr//g;
		($pp > 22) && next;
	
		($a[4] eq 'T' && $a[5] eq 'A') && next;
		($a[4] eq 'A' && $a[5] eq 'T') && next;
		($a[4] eq 'C' && $a[5] eq 'G') && next;
		($a[4] eq 'G' && $a[5] eq 'C') && next;
	

		print OUTF $a[0]."\t".($a[3]-1)."\t".$a[3]."\t".$a[1]."\n";
	}
	close(INF);
	close(OUTF);

	#- need to tabix it for bcftools below
	system("$bgzip -c ${resultsdir}/results/results_${sample}/temp/1kGP_GRCh38_nonAT_CG.bed > ${resultsdir}/results/results_${sample}/temp/1kGP_GRCh38_nonAT_CG.bed.gz");
	system("$tabix -p bed ${resultsdir}/results/results_${sample}/temp/1kGP_GRCh38_nonAT_CG.bed.gz");

}


sub _get_gvcf{

	#---------------------------------------
	#- get a gvcf file of the pgp sample that
	#-  contains only the 1kgp snps
	#-  and add rsIDs to the vcf file
	#---------------------------------------


    my $runstr="java -Xmx40g -Djava.io.tmpdir=${resultsdir}/results/results_${sample}/temp/tmpdir -jar $gatk";
    $runstr.=" -T HaplotypeCaller -R $ref_hs38";
    $runstr.=" -I $bam -nct $numThreads";
    $runstr.=" --emitRefConfidence GVCF -L ${resultsdir}/results/results_${sample}/temp/1kGP_GRCh38_nonAT_CG.bed";
    $runstr.=" -o ${resultsdir}/results/results_${sample}/temp/${sample}.g.vcf";

	#print "\n\nDEBUG NOT RUN: \n$runstr\n\n";
    `$runstr`;
		
	my $runstr2.="java -Xmx40g -Djava.io.tmpdir=${resultsdir}/results/results_${sample}/temp/tmpdir -jar $gatk";
	$runstr2.=" -T GenotypeGVCFs -R $ref_hs38";
	$runstr2.=" --variant ${resultsdir}/results/results_${sample}/temp/${sample}.g.vcf";
	$runstr2.=" -allSites -o ${resultsdir}/results/results_${sample}/temp/${sample}.genotypes.vcf";
	$runstr2.=" -stand_emit_conf 10 -stand_call_conf 30";
    `$runstr2`;
	
	system("cat ${resultsdir}/results/results_${sample}/temp/${sample}.genotypes.vcf | $bcftools annotate -a ${resultsdir}/results/results_${sample}/temp/1kGP_GRCh38_nonAT_CG.bed.gz -c CHROM,-,POS,ID | $bgzip -c > ${resultsdir}/results/results_${sample}/temp/${sample}.rsIDs.gvcf.gz");

}



sub _gvcf_to_plink() {

	#---------------------------------------
	#- convert gvcf to plink
	#---------------------------------------

	my $runstr1q="$plink";
	$runstr1q.=" --biallelic-only --vcf-require-gt";
	$runstr1q.=" --vcf ${resultsdir}/results/results_${sample}/temp/${sample}.rsIDs.gvcf.gz";
	$runstr1q.=" --out ${resultsdir}/results/results_${sample}/temp/$sample";
	$runstr1q.=' --allow-extra-chr --make-bed --double-id';	
	`$runstr1q`;

	#- the vcf has the SN (sample name) from the BAM, so when converts that to plink it gives me problems if there is a "_" in the name. So I'm changin it
	#	and I put the same $sample name as FIID and IID. Also, if the IID is not the same as $sample, the R script will not find the sameple and will not be plotted

	my $cmd="awk \'{\$1 = \"$sample\"; print}\' ${resultsdir}/results/results_${sample}/temp/${sample}.fam  > ${resultsdir}/results/results_${sample}/temp/${sample}.fam.mod";
	`$cmd`;
	my $cmd2="awk \'{\$2 = \"$sample\"; print}\' ${resultsdir}/results/results_${sample}/temp/${sample}.fam.mod  > ${resultsdir}/results/results_${sample}/temp/${sample}.fam";
	`$cmd2`;

}


sub _merge_pgp_1kGP {

	#---------------------------------------------
	#-- MERGE PGP SAMPLE WITH 1kGP GENOTYPES AND
	#---------------------------------------------


	#- merge the pgp and 1kGP files and get the SNPs that disagree between datasets
	
	my $runstr3="$plink";
	$runstr3.=" --bfile ${resultsdir}/results/results_${sample}/temp/${sample}";
	$runstr3.=" --bmerge ${initialAnc}";
	$runstr3.=" --out ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP_0";
	$runstr3.=" --geno 0 --allow-extra-chr";
	$runstr3.=" --make-bed";
	`$runstr3`;
	
	#- remove the discarded SNPs from the pgp and 1kGP file
	
	my $runstr1z="$plink";
	$runstr1z.=" --bfile ${initialAnc}";
	$runstr1z.=" --exclude ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP_0-merge.missnp";
	$runstr1z.=" --out ${resultsdir}/results/results_${sample}/temp/1kGP_2";
	$runstr1z.=" --make-bed --allow-extra-chr";
	`$runstr1z`;
	
	
	my $runstr1z2="$plink";
	$runstr1z2.=" --bfile ${resultsdir}/results/results_${sample}/temp/${sample}";
	$runstr1z2.=" --exclude ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP_0-merge.missnp";
	$runstr1z2.=" --out ${resultsdir}/results/results_${sample}/temp/${sample}_2";
	$runstr1z2.=" --make-bed --allow-extra-chr";
	`$runstr1z2`;
	
	
	#- merge the pgp and 1kGP files again
	
	my $runstr3a="$plink";
	$runstr3a.=" --bfile ${resultsdir}/results/results_${sample}/temp/${sample}_2";
	$runstr3a.=" --bmerge ${resultsdir}/results/results_${sample}/temp/1kGP_2";
	$runstr3a.=" --out ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP";
	$runstr3a.=" --geno 0 --allow-extra-chr";
	$runstr3a.=" --make-bed";
	`$runstr3a`;


}


sub _subset_unlinked_snps {

	#---------------------------------------------
	#-- CHOOSE A SUBSET OF FILTERED AND
	#--	INDEPENDENDT (UNLINKED) SNPS
	#---------------------------------------------


	#- select a subgroup of unlinked SNPs, by pruning those with r2 > 0.1 using 100-SNP windows
	#-	shifted at 5-SNP intervals. 
	
	my $runstr4="$plink";
	$runstr4.=" --bfile ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP";
	$runstr4.=" --out ${resultsdir}/results/results_${sample}/temp/snps_to_prune";
	$runstr4.=" --indep-pairwise 100 5 0.1";
	`$runstr4`;
	
	
	my $runstr5="$plink";
	$runstr5.=" --bfile ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP";
	$runstr5.=" --out ${resultsdir}/results/results_${sample}/temp/${sample}_1kGP_pruned";
	$runstr5.=" --extract ${resultsdir}/results/results_${sample}/temp/snps_to_prune.prune.in";
	$runstr5.=" --make-bed --mind 0.1";
	`$runstr5`;

}


