#!/usr/bin/env python3

import os
import subprocess
from pathlib import Path

import fire


def get_list_of_snps_from_1kGP(initialBIM, output_dir, hasCHRflag=0, bgzip='bgzip', tabix='tabix'):
    """
    Get list of SNPs from 1kGP

    Input files:
        - initialBIM = f"{dir}reference/1kGP_GRCh38_exome.bim"
    Output files:
        - f"{dir}temp/1kGP_GRCh38_nonAT_CG.bed.gz"
        - f"{dir}temp/1kGP_GRCh38_nonAT_CG.bed.gz.tbi"

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir/f"1kGP_GRCh38_nonAT_CG.bed"
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
    os.system(f"{bgzip} -c {output_dir}/1kGP_GRCh38_nonAT_CG.bed > {output_dir}/1kGP_GRCh38_nonAT_CG.bed.gz")
    os.system(f"{tabix} -p bed {output_dir}/1kGP_GRCh38_nonAT_CG.bed.gz")


def get_gvcf(vcf_or_bam, output_dir, ref_hs38, ref_1kGP, name='ancestry',
             gatk='software/GenomeAnalysisTK.jar', numThreads=4, bcftools='bcftools', bgzip='bgzip'):
    """
    Convert BAM to a gvcf that contains only the 1kgp snps filtered above

    Input files:
        - f"{dir}bam/clean.bam"

        - {ref_hs38} reference/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa
        - {ref_1kGP} {tmp_dir}/1kGP_GRCh38_nonAT_CG.bed

    Output files:
        - f"{dir}temp/g.vcf"
        - f"{dir}temp/genotypes.vcf"
        - f"{dir}temp/rsIDs.gvcf.gz"

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    ref_1kGP = Path(ref_1kGP)
    ref_1kGP_gz = f"{ref_1kGP}.gz"

    if str(vcf_or_bam).lower().endswith('.bam'):
        # Define the command to run GATK HaplotypeCaller
        output_g_path = f"{output_dir}/{name}.g.vcf"
        if not Path(output_g_path).exists():
            runstr = f"java -Xmx40g -Djava.io.tmpdir={output_dir}/tmpdir -jar {gatk}"
            runstr += f" -T HaplotypeCaller -R {ref_hs38}"
            runstr += f" -I {vcf_or_bam} -nct {numThreads}"
            runstr += f" --emitRefConfidence GVCF -L {ref_1kGP}"
            runstr += f" -o {output_g_path}"
            os.system(runstr)

        # Define the command to run GATK GenotypeGVCFs
        output_genotypes_path = f"{output_dir}/{name}.genotypes.vcf"
        if not Path(output_genotypes_path).exists():
            # return
            runstr2 = f"java -Xmx40g -Djava.io.tmpdir={output_dir}/tmpdir -jar {gatk}"
            runstr2 += f" -T GenotypeGVCFs -R {ref_hs38}"
            runstr2 += f" --variant {output_g_path}"
            runstr2 += f" -allSites -o {output_genotypes_path}"
            runstr2 += " -stand_emit_conf 10 -stand_call_conf 30"
            os.system(runstr2)
    else:
        output_genotypes_path = vcf_or_bam

    # Add rsIDs to the VCF file using bcftools and bgzip
    output_path = f"{output_dir}/{name}.rsIDs.gvcf.gz"
    if not Path(output_path).exists():
        cmd = f"cat {output_genotypes_path} | {bcftools} annotate -a {ref_1kGP_gz} -c CHROM,-,POS,ID | {bgzip} -c > {output_path}"
        ret = subprocess.run(cmd, shell=True)


def gvcf_to_plink(gvcf_path, output_dir, output_name='sample', plink='software/plink'):
    """
    Convert gvcf to plink files (bed, bim, fam)

    Input files:
        - {output_dir}/{sample}.ancestry.rsIDs.gvcf.gz

    Output files:
        - {output_dir}/{sample}.bed
        - {output_dir}/{sample}.bim
        - {output_dir}/{sample}.fam
        - {output_dir}/{sample}.fam.mod
        - {output_dir}/{sample}.log (update)
        - {output_dir}/{sample}.nosex

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = f"{output_dir}/{output_name}.fam"
    if Path(output_path).exists():
        return

    # convert gvcf to plink (generate bed, bim, fam)
    runstr1q = f"{plink} --biallelic-only --vcf-require-gt"
    runstr1q += f" --vcf {gvcf_path}"
    runstr1q += f" --out {output_dir}/{output_name}"
    runstr1q += " --allow-extra-chr --make-bed --double-id"
    subprocess.run(runstr1q, shell=True)

    # the vcf has the SN (sample name) from the BAM, so when converts that to plink it gives me problems if there is a "_" in the name.
    # So I'm changin it and I put the same $sample name as FIID and IID.
    # Also, if the IID is not the same as $sample, the R script will not find the sameple and will not be plotted
    cmd = f"awk '{{ $1 = \"{output_name}\"; print }}' {output_dir}/{output_name}.fam > {output_dir}/{output_name}.fam.mod"
    subprocess.run(cmd, shell=True)
    cmd2 = f"awk '{{ $2 = \"{output_name}\"; print }}' {output_dir}/{output_name}.fam.mod > {output_dir}/{output_name}.fam"
    subprocess.run(cmd2, shell=True)


def merge_pgp_1kGP(plink_result_dir, output_dir, ref_initialAnc='reference/1kGP_GRCh38_exome', plink='software/plink'):
    """
    Merge the pgp and 1kGP files and get the SNPs that disagree between datasets
    - Four groups of files are generated:

    Input files:
        - {plink_result_dir}/{sample}.bed
        - {plink_result_dir}/{sample}.bim
        - {plink_result_dir}/{sample}.fam
        - {plink_result_dir}/{sample}.fam.mod
        - initialAnc f"{dir}reference/1kGP_GRCh38_exome"

    Output files:
        - {resultsdir}/{sample}_1kGP_0 (-merge.missnp, -merge.fam, .log)
        - {resultsdir}/1kGP_2 (.bed, .bim, .fam, .log, .nosex)
        - {resultsdir}/{sample}_2 (.bed, .bim, .fam, .log, .nosex)
        - {resultsdir}/{sample}_1kGP (.bed, .bim, .fam, .log, .nosex)

    """
    # reusltsdir = f"{resultsdir}/temp"
    # output_path = f"{reusltsdir}/{sample}.fam"
    # if Path(output_path).exists():
    #     return
    # Merge the pgp and 1kGP files and get the SNPs that disagree between datasets
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    output_path = f"{output_dir}/1kGP_0-merge.fam"
    if not Path(output_path).exists():
        runstr3 = f"{plink} --bfile {plink_result_dir} --bmerge {ref_initialAnc} --out {output_dir}/1kGP_0 --geno 0 --allow-extra-chr --make-bed"
        os.system(runstr3)
    
    # Remove discarded SNPs from the pgp and 1kGP file
    output_path = f"{output_dir}/1kGP_2.bed"
    if not Path(output_path).exists():
        runstr1z = f"{plink} --bfile {ref_initialAnc} --exclude {output_dir}/1kGP_0-merge.missnp --out {output_dir}/1kGP_2 --make-bed --allow-extra-chr"
        os.system(runstr1z)

    output_path = f"{output_dir}/2.bed"
    if not Path(output_path).exists():
        runstr1z2 = f"{plink} --bfile {plink_result_dir} --exclude {output_dir}/1kGP_0-merge.missnp --out {output_dir}/2 --make-bed --allow-extra-chr"
        os.system(runstr1z2)
    
    # Merge the pgp and 1kGP files again
    output_path = f"{output_dir}/1kGP.bed"
    if not Path(output_path).exists():
        runstr3a = f"{plink} --bfile {output_dir}/2 --bmerge {output_dir}/1kGP_2 --out {output_dir}/1kGP --geno 0 --allow-extra-chr --make-bed"
        os.system(runstr3a)


def subset_unlinked_snps(plink_1kgp_path, output_dir, plink='software/plink'):
    """
    Choose a subset of filtered and independent (unlinked) SNPs

    Input files:
        - {output_dir}/{sample}_1kGP* (.bed, .bim, .fam, .log, .nosex)

    Output files:
        - {output_dir}/snps_to_prune (.log, .prune.in, .prune.out, .nosex)
        - {output_dir}/{sample}_1kGP_pruned (.bed, .bim, .fam, .log, .nosex)

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # output_dir = f"{resultsdir}/temp"

    # Choose a subset of filtered and independent (unlinked) SNPs

    # Select a subgroup of unlinked SNPs, by pruning those with r2 > 0.1 using 100-SNP windows shifted at 5-SNP intervals.
    output_path = f"{output_dir}/snps_to_prune.prune.in"
    if not Path(output_path).exists():
        runstr4 = f"{plink} --bfile {plink_1kgp_path} --out {output_dir}/snps_to_prune --indep-pairwise 100 5 0.1"
        os.system(runstr4)

    # Remove the pruned SNPs from the dataset and set the missingness threshold to 0.1
    output_path = f"{output_dir}/1kGP_pruned.bed"
    if not Path(output_path).exists():
        runstr5 = f"{plink} --bfile {plink_1kgp_path} --out {output_dir}/1kGP_pruned --extract {output_dir}/snps_to_prune.prune.in --make-bed --mind 0.1"
        os.system(runstr5)

def from_1kgp_pruned_to_pca(plink_1kgp_pruned, output_dir, plink='software/plink'):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    output_path = f"{output_dir}/1kGP_pruned_pca_20.eigenvec"
    if not Path(output_path).exists():
        runstr6 = f"{plink}"
        runstr6 += f" --bfile {plink_1kgp_pruned}"
        runstr6 += f" --out {output_dir}/1kGP_pruned_pca_20"
        runstr6 += " --pca"
        subprocess.run(runstr6, shell=True)

def generate_ancestry_pipeline(bam_or_vcf, output_dir, numThreads=4, sample='sample',
                               dir_reference='reference', dir_software='software', workdir="",
                               ):
    """
    Generate ancestry from BAM file

    Steps:
        1. preprocessed the 1kGP data
        2. Convert BAM to VCF
        3. Convert VCF to PLINK
        4. Merge the pgp and 1kGP files and get the SNPs that disagree between datasets
        5. Choose a subset of filtered and independent (unlinked) SNPs
        6. Run ADMIXTURE


    Parameters
    ----------
    bam
    output_dir
    numThreads
    dir_reference
    dir_software
    workdir

    Returns
    -------

    """

    dir = ""
    if 'SINGULARITY_NAME' in os.environ:
        dir = "/GenomeChronicler/"

    # resultsdir = dir
    plink = f"{dir}software/plink"
    gatk = f"{dir}software/GenomeAnalysisTK.jar"
    ref_hs38 = f"{dir}reference/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa"
    initialBIM = f"{dir}reference/1kGP_GRCh38_exome.bim"
    initialAnc = f"{dir}reference/1kGP_GRCh38_exome"
    bcftools = "bcftools"
    bgzip = "bgzip"
    tabix = "tabix"
    hasCHRflag = 0

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    # Also check for the need to samtools index this stuff
    # java -jar software/picard.jar CreateSequenceDictionary R=software/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa O=software/GRCh38_full_analysis_set_plus_decoy_hla_noChr.dict

    # ------------------------------------------
    # -- I had already preprocessed the 1kGP data to reduce it to a PLINK file with
    # -- only SNPs from exome (914577 total)
    # -- But from these, here I remove all the A/T or C/G SNPs to make sure that there
    # -- will not be artifacts due to reporting SNPS i different stransd ib both datasets
    # -- I also remove sexual chromosomes
    # ------------------------------------------
    plink_1kgp_dir = output_dir/'1kGP_GRCh38_nonAT_CG'
    get_list_of_snps_from_1kGP(initialBIM, plink_1kgp_dir, hasCHRflag=hasCHRflag, bgzip=bgzip, tabix=tabix)

    # ------------------------------------------
    # -- Convert BAM to a gvcf that contains only
    # - the 1kgp snps filtered above
    # ------------------------------------------
    ref_1kGP = plink_1kgp_dir/'1kGP_GRCh38_nonAT_CG.bed'
    vcf_dir = output_dir/'vcf_ancestry'
    get_gvcf(bam_or_vcf, vcf_dir, ref_hs38, ref_1kGP, 'ancestry', gatk, numThreads, bcftools, bgzip)

    # ------------------------------------------
    # -- Convert GVCF data to PLINK format
    # ------------------------------------------
    gvcf_path = vcf_dir/f"ancestry.rsIDs.gvcf.gz"
    gvcf_to_plink(gvcf_path, vcf_dir, sample, plink)

    # ------------------------------------------
    # -- Merge gvcf with 1kGP data. I need to merge the files twice, the first one will fail if there are SNPs with
    # -- different alleles in both datasets. The second merge will exclude these SNPs
    # ------------------------------------------
    plink_result_dir = vcf_dir/sample
    merge_pgp_1kgp_dir = output_dir/'merged_pgp_1kGP'
    merge_pgp_1kGP(plink_result_dir, merge_pgp_1kgp_dir, initialAnc, plink)

    # ------------------------------------------
    # -- Reduce the merged dataset further by removing SPs in LD
    # ------------------------------------------
    merge_pgp_1kgp_path = merge_pgp_1kgp_dir/f"1kGP"
    subset_unlinked_snps_dir = output_dir/'subset_unlinked_snps'
    subset_unlinked_snps(merge_pgp_1kgp_path, subset_unlinked_snps_dir, plink)

    # ------------------------------------------
    # -- Perform PCA
    # ------------------------------------------
    plink_result_dir = subset_unlinked_snps_dir/f"1kGP_pruned"
    pca_dir = output_dir/'pca_ancestry'
    from_1kgp_pruned_to_pca(plink_result_dir, pca_dir, plink)


funcs = {
    'get_list_of_snps_from_1kGP': get_list_of_snps_from_1kGP,
    'get_gvcf': get_gvcf,
    'gvcf_to_plink': gvcf_to_plink,
    'merge_pgp_1kGP': merge_pgp_1kGP,
    'subset_unlinked_snps': subset_unlinked_snps,
    'from_1kgp_pruned_to_pca': from_1kgp_pruned_to_pca,
    'generate_ancestry_pipeline': generate_ancestry_pipeline,
}


class Tools(object):
    def __init__(self):
        super(Tools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(Tools, k, staticmethod(v))
    fire.Fire(Tools)
