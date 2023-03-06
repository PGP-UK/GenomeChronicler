#!/usr/bin/env python3

import os
import sys
import subprocess
from os.path import basename

from pathlib import Path
import fire

# get list of SNPs from 1kGP
def _get_list_of_snps_from_1kGP(initialBIM, sample, resultsdir, hasCHRflag, bgzip, tabix):
    output_path = f"{resultsdir}/results/results_{sample}/temp/1kGP_GRCh38_nonAT_CG.bed"
    if Path(output_path).exists():
        return
    with open(initialBIM, "r") as inf, open(output_path, "w") as outf:
        outf.write("#CHROM\tST0\tPOS\tID\n")
        for line in inf:
            a = line.strip().split("\t")
            if hasCHRflag == 1 and not a[0].startswith("chr"):
                a[0] = "chr" + a[0]
            pp = a[0]
            pp = pp.replace('chr','')
            if int(pp) > 22:
                continue
            if a[4] == 'T' and a[5] == 'A':
                continue
            if a[4] == 'A' and a[5] == 'T':
                continue
            if a[4] == 'C' and a[5] == 'G':
                continue
            if a[4] == 'G' and a[5] == 'C':
                continue
            outf.write(f"{a[0]}\t{int(a[3])-1}\t{a[3]}\t{a[1]}\n")
    os.system(f"{bgzip} -c {resultsdir}/results/results_{sample}/temp/1kGP_GRCh38_nonAT_CG.bed > {resultsdir}/results/results_{sample}/temp/1kGP_GRCh38_nonAT_CG.bed.gz")
    os.system(f"{tabix} -p bed {resultsdir}/results/results_{sample}/temp/1kGP_GRCh38_nonAT_CG.bed.gz")


# convert BAM to a gvcf that contains only the 1kgp snps filtered above
def _get_gvcf(vcf, resultsdir, sample, bcftools, bgzip):
    tmp_dir = f"{resultsdir}/results/results_{sample}/temp"
    # output_path = f"{tmp_dir}/{sample}.rsIDs.gvcf.gz"
    # if Path(output_path).exists():
    #     return
    # runstr = f"java -Xmx40g -Djava.io.tmpdir={tmp_dir}/tmpdir -jar {gatk}"
    # runstr += f" -T HaplotypeCaller -R {ref_hs38}"
    # runstr += f" -I {bam} -nct {numThreads}"
    # runstr += f" --emitRefConfidence GVCF -L {tmp_dir}/1kGP_GRCh38_nonAT_CG.bed"
    # runstr += f" -o {tmp_dir}/{sample}.g.vcf"
    # os.system(runstr)
    # runstr2 = f"java -Xmx40g -Djava.io.tmpdir={tmp_dir}/tmpdir -jar {gatk}"
    # runstr2 += f" -T GenotypeGVCFs -R {ref_hs38}"
    # runstr2 += f" --variant {tmp_dir}/{sample}.g.vcf"
    # runstr2 += f" -allSites -o {tmp_dir}/{sample}.genotypes.vcf"
    # runstr2 += " -stand_emit_conf 10 -stand_call_conf 30"
    # os.system(runstr2)
    os.system(f"cat {vcf} | {bcftools} annotate -a {tmp_dir}/1kGP_GRCh38_nonAT_CG.bed.gz -c CHROM,-,POS,ID | {bgzip} -c > {tmp_dir}/{sample}.rsIDs.gvcf.gz")


def _gvcf_to_plink(plink, resultsdir, sample):
    tmp_dir = f"{resultsdir}/results/results_{sample}/temp"
    output_path = f"{tmp_dir}/{sample}.fam"
    if Path(output_path).exists():
        return
    # convert gvcf to plink
    runstr1q = f"{plink} --biallelic-only --vcf-require-gt"
    runstr1q += f" --vcf {tmp_dir}/{sample}.rsIDs.gvcf.gz"
    runstr1q += f" --out {tmp_dir}/{sample}"
    runstr1q += " --allow-extra-chr --make-bed --double-id"
    subprocess.run(runstr1q, shell=True)

    # the vcf has the SN (sample name) from the BAM, so when converts that to plink it gives me problems if there is a "_" in the name.
    # So I'm changin it and I put the same $sample name as FIID and IID.
    # Also, if the IID is not the same as $sample, the R script will not find the sameple and will not be plotted
    cmd = f"awk '{{ $1 = \"{sample}\"; print }}' {tmp_dir}/{sample}.fam > {tmp_dir}/{sample}.fam.mod"
    subprocess.run(cmd, shell=True)
    cmd2 = f"awk '{{ $2 = \"{sample}\"; print }}' {tmp_dir}/{sample}.fam.mod > {tmp_dir}/{sample}.fam"
    subprocess.run(cmd2, shell=True)

# def _gvcf_to_plink(plink, resultsdir, sample):
#     tmp_dir = f"{resultsdir}/results/results_{sample}/temp"
#     # Convert gvcf to plink
#     runstr1q = f"{plink} --biallelic-only --vcf-require-gt --vcf {tmp_dir}/{sample}.rsIDs.gvcf.gz --out {tmp_dir}/{sample} --allow-extra-chr --make-bed --double-id"
#     os.system(runstr1q)
    
#     # Change sample names to avoid errors and inconsistencies
#     cmd = f"awk '{{\$1 = \"{sample}\"; print}}' {tmp_dir}/{sample}.fam  > {tmp_dir}/{sample}.fam.mod"
#     os.system(cmd)
#     cmd2 = f"awk '{{\$2 = \"{sample}\"; print}}' {tmp_dir}/{sample}.fam.mod  > {tmp_dir}/{sample}.fam"
#     os.system(cmd2)


def _merge_pgp_1kGP(plink, resultsdir, sample, initialAnc):
    tmp_dir = f"{resultsdir}/results/results_{sample}/temp"
    # output_path = f"{tmp_dir}/{sample}.fam"
    # if Path(output_path).exists():
    #     return
    # Merge the pgp and 1kGP files and get the SNPs that disagree between datasets
    runstr3 = f"{plink} --bfile {tmp_dir}/{sample} --bmerge {initialAnc} --out {tmp_dir}/{sample}_1kGP_0 --geno 0 --allow-extra-chr --make-bed"
    os.system(runstr3)
    
    # Remove discarded SNPs from the pgp and 1kGP file
    runstr1z = f"{plink} --bfile {initialAnc} --exclude {tmp_dir}/{sample}_1kGP_0-merge.missnp --out {tmp_dir}/1kGP_2 --make-bed --allow-extra-chr"
    os.system(runstr1z)
    runstr1z2 = f"{plink} --bfile {tmp_dir}/{sample} --exclude {tmp_dir}/{sample}_1kGP_0-merge.missnp --out {tmp_dir}/{sample}_2 --make-bed --allow-extra-chr"
    os.system(runstr1z2)
    
    # Merge the pgp and 1kGP files again
    runstr3a = f"{plink} --bfile {tmp_dir}/{sample}_2 --bmerge {tmp_dir}/1kGP_2 --out {tmp_dir}/{sample}_1kGP --geno 0 --allow-extra-chr --make-bed"
    os.system(runstr3a)


def _subset_unlinked_snps(plink, resultsdir, sample):
    tmp_dir = f"{resultsdir}/results/results_{sample}/temp"
    output_path = f"{tmp_dir}/{sample}_1kGP_pruned.bed"
    if Path(output_path).exists():
        return
    
    # Choose a subset of filtered and independent (unlinked) SNPs
    
    # Select a subgroup of unlinked SNPs, by pruning those with r2 > 0.1 using 100-SNP windows shifted at 5-SNP intervals.
    runstr4 = f"{plink} --bfile {tmp_dir}/{sample}_1kGP --out {tmp_dir}/snps_to_prune --indep-pairwise 100 5 0.1"
    os.system(runstr4)
    
    # Remove the pruned SNPs from the dataset and set the missingness threshold to 0.1
    runstr5 = f"{plink} --bfile {tmp_dir}/{sample}_1kGP --out {tmp_dir}/{sample}_1kGP_pruned --extract {tmp_dir}/snps_to_prune.prune.in --make-bed --mind 0.1"
    os.system(runstr5)


def ancestry_generator_from_VCF(vcf, resultsdir=""):

    dir = ""
    if 'SINGULARITY_NAME' in os.environ:
        dir = "/GenomeChronicler/"

    # resultsdir = dir
    plink = f"{dir}software/plink"
    # gatk = f"{dir}software/GenomeAnalysisTK.jar"
    # ref_hs38 = f"{dir}reference/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa"
    initialBIM = f"{dir}reference/1kGP_GRCh38_exome.bim"
    initialAnc = f"{dir}reference/1kGP_GRCh38_exome"
    bcftools = "bcftools"
    # samtools = "samtools"
    bgzip = "bgzip"
    tabix = "tabix"

    hasCHRflag = 0
    numThreads = 4

    # if not sys.argv:
    #     print("Please provide a single base recalibrated VCF file to this script")
    #     sys.exit(1)

    # VCF = sys.argv[0]
    sample = basename(vcf)
    sample = sample.rsplit(".", 1)[0]
    sample = sample.replace(".recal", "")
    # sample = sample.replace(".bam.clean", "")

    # if len(sys.argv) > 1:
    #     resultsdir = sys.argv[1]

    # if len(sys.argv) > 2:
    #     numThreads = sys.argv[2]

    os.system(f"mkdir -p {resultsdir}/results/results_{sample}/temp")

    # Also check for the need to samtools index this stuff
    # java -jar software/picard.jar CreateSequenceDictionary R=software/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa O=software/GRCh38_full_analysis_set_plus_decoy_hla_noChr.dict

    # ------------------------------------------
    # -- I had already preprocessed the 1kGP data to reduce it to a PLINK file with
    # -- only SNPs from exome (914577 total)
    # -- But from these, here I remove all the A/T or C/G SNPs to make sure that there
    # -- will not be artifacts due to reporting SNPS i different stransd ib both datasets
    # -- I also remove sexual chromosomes
    # ------------------------------------------

    # _get_list_of_snps_from_1kGP()
    _get_list_of_snps_from_1kGP(initialBIM, sample, resultsdir, hasCHRflag, bgzip, tabix)

    # ------------------------------------------
    # -- Convert BAM to a gvcf that contains only
    # - the 1kgp snps filtered above
    # ------------------------------------------

    # _get_gvcf()
    # _get_gvcf(vcf, resultsdir, sample, gatk, ref_hs38, numThreads, bcftools, bgzip)
    _get_gvcf(vcf, resultsdir, sample, bcftools, bgzip)

    # ------------------------------------------
    # -- Convert GVCF data to PLINK format
    # ------------------------------------------

    # _gvcf_to_plink()
    _gvcf_to_plink(plink, resultsdir, sample)

    # ------------------------------------------
    # -- Merge gvcf with 1kGP data. I need to merge the files twice, the first one will fail if there are SNPs with
    # -- different alleles in both datasets. The second merge will exclude these SNPs
    # ------------------------------------------

    # _merge_pgp_1kGP()
    _merge_pgp_1kGP(plink, resultsdir, sample, initialAnc)
    # print('test')
    # exit()

    # ------------------------------------------
    # -- Reduce the merged dataset further by removing SPs in LD
    # ------------------------------------------

    # _subset_unlinked_snps()
    _subset_unlinked_snps(plink, resultsdir, sample)
    # ------------------------------------------
    # -- Perform PCA
    # ------------------------------------------
    tmp_dir = f"{resultsdir}/results/results_{sample}/temp"
    output_path = f"{tmp_dir}/{sample}_1kGP_pruned_pca_20.eigenvec"
    if Path(output_path).exists():
        return
    runstr6 = f"{plink}"
    runstr6 += f" --bfile {tmp_dir}/{sample}_1kGP_pruned"
    runstr6 += f" --out {tmp_dir}/{sample}_1kGP_pruned_pca_20"
    runstr6 += " --pca"
    subprocess.run(runstr6, shell=True)


if __name__ == '__main__':
    fire.Fire(ancestry_generator_from_VCF)