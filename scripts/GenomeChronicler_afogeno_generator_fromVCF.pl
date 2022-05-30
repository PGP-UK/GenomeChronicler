#!/usr/bin/env perl

use warnings;
use strict;
use File::Basename;


my $dir="/GenomeChronicler/";
my $resultsdir = $dir;

my $inputBED = "${dir}/reference/snps.19-114.unique.nochr.bed";
my $gatk="${dir}/software/GenomeAnalysisTK.jar";
#my $ref_hs37="${dir}/software/human_g1k_v37_decoy.fasta";
my $ref_hs38="${dir}/reference/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa";
# my $bcftools="${dir}/software/bcftools";
# my $samtools="${dir}/software/samtools";
# my $bgzip="${dir}software/bgzip";
# my $tabix="${dir}software/tabix";
my $bcftools="bcftools";
my $samtools="samtools";
my $bgzip="bgzip";
my $tabix="tabix";


my $numThreads = 4;


#Needs tabix and bgzip here in the tools too; 

if(!@ARGV) {
	die "Please provide a single base recalibrated gVCF file to this script and optional output directory\n";
}

my $gVCF=$ARGV[0];
(my $sample= basename($gVCF)) =~ s/\.[^.]+$//;
$sample =~ s/\.recal//g;
$sample =~ s/\.gVCF\.clean//gi;

$resultsdir = $ARGV[1] if(defined($ARGV[1]));

$numThreads = $ARGV[2] if(defined($ARGV[2]));

system("mkdir -p ${resultsdir}/results/results_${sample}/temp");

#die "Try to get all iID coordinates from all 23andme chip versions - ok for 37, but can we get it for 38 now?\n";


#Check if gVCF exists
if(!-e $gVCF) {
	die "The gVCF file specified does not exist, please double-check the path and try again [$gVCF]\n";
}


#Check if all BED files are in place;
if(!-e $inputBED) {
	die "The BED file specified does not exist, please double-check the path and try again [$inputBED]\n";
}


if(!-e ${inputBED}.".gz") {

	#Sort and compress BED
	system("sort -k 1,1 -k2,2n $inputBED > tempAAA");
	system("mv tempAAA $inputBED");
	system("$bgzip -c ${inputBED} > ${inputBED}.gz");

	#Force tabix re-index
	system("$tabix -p bed ${inputBED}.gz");
}

if(!-e ${inputBED}.".gz.tbi") {
	#tabix index
	system("$tabix -p bed ${inputBED}.gz");
}


sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}



#---------------------------------------
#- get a vcf file of the pgp sample that
#-  contains only the above snps
#---------------------------------------

#Might be useful after JAVA : -Xmx40g -Djava.io.tmpdir=${resultsdir}/results/results_${sample}/temp/tmpdir 


# -o ${resultsdir}/results/results_${sample}/temp/${sample}.genotypes.vcf";


#---------------------------------------
#- add rsIDs to the vcf file
#---------------------------------------

system("cat ${resultsdir}/results/results_${sample}/temp/${sample}.genotypes.vcf | $bcftools annotate -a ${inputBED}.gz -c CHROM,-,POS,ID | $bgzip -c > ${resultsdir}/results/results_${sample}/temp/${sample}.genotypes.rsIDs.vcf.gz");


##################### Generate AFO Geno 38 file


#Now that we have a bed file, read it and keep it somewhere safe, it is stupid to read it over and over again...
my %genData;
my %counterTemplate;

open(IN, "gzip -dcf $inputBED |") or die "Could not open input file: $!\n";

while (my $line = <IN>) {
    chomp($line);
    next if($line eq "");
    my @data = split("\t", $line);
    
    #$data[0] = "chr".$data[0] if($data[0] !~ m/^chr/);
    #die "@data\n";
    
    $genData{$data[0]}{$data[2]} = \@data;
    $counterTemplate{$data[0]}{$data[2]} = 0;
 }

close IN;

#die Dumper(\%genData);

my $VCFfilename  = "${resultsdir}/results/results_${sample}/temp/${sample}.genotypes.rsIDs.vcf.gz";
my %debugCounter = %counterTemplate;

open(IN, "gzip -dcf $VCFfilename | grep -v \"0/0:0:0:0:0,0,0\" | grep -v LowQual | cut -f 1,2,4,5,10 | uniq |") or die "Could not open input file: $!\n";
open(OUT, ">${resultsdir}/results/results_${sample}/temp/${sample}.afogeno38.txt") or die "Could not open output file: $!\n";

while (my $line = <IN>) {
	chomp($line);

    #print STDERR "SUPERDEBUG: $line\n";

	next if($line eq "");
	next if($line =~ m/^#/);
        
    $line =~ s/,*\<NON_REF\>//g;
        
	my @data = split("\t", $line);

	my @tmp = split(":", $data[4]);
	$data[4] = $tmp[0];

	#note: do something about the ./. here?
	my ($genotype,$extra) = (["NA","NA"],"NA");

	if($data[4] ne "./.") {
		($genotype,$extra) =  buildGenotype($data[2],$data[3],$data[4]);
	}
	else { #No reads there
        $debugCounter{$data[0]}{$data[1]} += 1;
		next;
	}

	push(@data,@{$genotype});
	push(@data,$extra);

    #die Dumper(\@data);
    if(exists($genData{$data[0]}{$data[1]})) {
        $debugCounter{$data[0]}{$data[1]} += 1;
        my $d2 = $genData{$data[0]}{$data[1]};
        
        local $" ="\t";
        print OUT "@{$d2}\t@data\n";
    }
    
}

close IN;
close OUT;


#Ensuring VCF gets stored for posterity like requested by reviewer 4. 
system ("mv $VCFfilename ${resultsdir}/results/results_${sample}/${sample}.genotypingVCF.vcf.gz");

# print STDERR "Finished processing file, now for the debug routines...\n";

### Commented to save space when running on AWS - re-enable for debug
# for my $key1 (keys(%debugCounter)) {
#     for my $key2(keys(%{$debugCounter{$key1}})) {
            
#         #In general check if debug is exactly 1, any deviations should be reported. (0 = not in gvcf; 2 = repeated in gvcf and a sign of trouble);
#         print STDERR "+++ ERROR\tLikely LowQual counts\t[$debugCounter{$key1}{$key2}]\t--\t$key1\t$key2\t$genData{$key1}{$key2}[3]\t+++\n" if($debugCounter{$key1}{$key2} != 1);
#     }
# }








sub buildGenotype {
    my @alleles = ($_[0],split(",",$_[1]));
    my @extras = ("");
    my $gen = $_[2];
    
    
    my $l1 = length($alleles[0]);
    my $extra = "";
    

    for (my $i = 1; $i < @alleles; $i++ ) {
        
        $extras[$i] = "";
        
        my $l2 = length($alleles[$i]);
        if($l1 > 1 or $l2 > 1) {
            
            ##Debug statement
            #print STDERR "There is an indel here [@alleles] [$gen]\n";
            
            
            if($l2 > $l1) {
                $alleles[0] = "-";
                $alleles[$i] =~ s/^.//s ;
                #die "There is an insertion here [@alleles] [$gen]\n";
                $extras[$i] = "I";
            }
            else {
                $alleles[0] =~ s/^.//s ;
                $alleles[$i] = "-";
                #die "Implement deletions here please  [@alleles] [$gen]\n";
                $extras[$i] = "D";
            }
        }
    }

    
    #    die "Genotype debug [$gen]\n" if($gen ne "0/1" and $gen ne "1/1"); #This is to account for odd genotypes. None have been encountered so far but they exist... #Now solved in this version I think
    
    my @res = map{$alleles[$_]} split('/',$gen);
    $extra = join(";",map{$extras[$_]} split('/',$gen));
    $extra =~ s/^;//;
    $extra =~ s/;$//;
    
    #print STDERR "@alleles\t$gen\t@res\n";
        
    #die "Just come back here to buildGenotype and check this case is working as it should\n=====\nINPUT:@_ -- @alleles -- $gen == @res --- $extra\n=====\n" if(length($extra) > 1);
    
    #die "Just come back here to buildGenotype and check this case is working as it should\n=====\nINPUT:@_ -- [[@alleles]] -- $gen == @res --- $extra\n=====\n" if($_[2] =~ m/0/);
    
    
    print STDERR "\n=====\nINPUT:@_ -- @alleles -- $gen == @res\n=====\n" if($_[2] =~ m/2/);
    
    return \@res,$extra;
}



