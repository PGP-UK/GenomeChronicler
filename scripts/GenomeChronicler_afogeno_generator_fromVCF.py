#!/usr/bin/env python3
import subprocess
import os
import sys
import gzip
from os.path import basename
from tqdm import tqdm
from pathlib import Path

import fire


def build_genotype(a1, a2, gen):
    alleles = [a1] + a2.split(",")
    extras = [""] * len(alleles)
    l1 = len(alleles[0])
    extra = ""
    for i in range(1, len(alleles)):
        extras[i] = ""
        l2 = len(alleles[i])
        if l1 > 1 or l2 > 1:
            if l2 > l1:
                alleles[0] = "-"
                alleles[i] = alleles[i][1:]
                extras[i] = "I"
            else:
                alleles[0] = alleles[0][1:]
                alleles[i] = "-"
                extras[i] = "D"
    res = [alleles[int(i)] for i in gen.split("/")]
    extra = ";".join([extras[int(i)] for i in gen.split("/")]).strip(";")
    return res, extra


def uniq(l):
    return list(set(l))


def afogeno_generator_from_VCF(input_VCF, sample, resultsdir="",
                               dir_reference='reference', dir_software='software', workdir="", ):
    if 'SINGULARITY_NAME' in os.environ:
        workdir = "/GenomeChronicler/"

    ref_bed = workdir + f"{dir_reference}/snps.19-114.unique.nochr.bed"
    bcftools = "bcftools"
    bgzip = workdir + f"{dir_software}/bgzip"
    tabix = workdir + f"{dir_software}/tabix"

    # sample = basename(input_VCF).rsplit(".", 1)[0]

    subprocess.run(["mkdir", "-p", f"{resultsdir}/temp"])

    if not os.path.exists(input_VCF):
        sys.exit(f"The VCF file specified does not exist, please double-check the path and try again [{input_VCF}]")

    if not os.path.exists(ref_bed):
        sys.exit(f"The BED file specified does not exist, please double-check the path and try again [{ref_bed}]")

    if not os.path.exists(ref_bed + ".gz"):
        # Sort and compress BED
        subprocess.run(["sort", "-k", "1,1", "-k2,2n", ref_bed, ">", "tempAAA"])
        subprocess.run(["mv", "tempAAA", ref_bed])
        subprocess.run([bgzip, "-c", f"{ref_bed}" + " > " + f"{ref_bed}.gz"])

        # Force tabix re-index
        subprocess.run([tabix, "-p", "bed", f"{ref_bed}.gz"])

    if not os.path.exists(f"{ref_bed}.gz.tbi"):
        # tabix index
        subprocess.run([tabix, "-p", "bed", f"{ref_bed}.gz"])

    # Add rsIDs to the VCF file using bcftools and bgzip
    # os.system(f"cat {resultsdir}/temp/{sample}.genotypes.vcf | {bcftools} annotate -a {ref_bed}.gz -c CHROM,-,POS,ID | {bgzip} -c > {resultsdir}/temp/{sample}.genotypes.rsIDs.vcf.gz")
    # vcf_path = f"{resultsdir}/temp/{sample}.genotypes.vcf"

    vcf_path = input_VCF
    bed_path = f"{ref_bed}.gz"
    out_path = f"{resultsdir}/temp/{sample}.genotypes.rsIDs.vcf.gz"

    # subprocess.run(f"cat {vcf_path} | {bcftools} annotate -a {bed_path} -c CHROM,-,POS,ID | {bgzip} -c > {out_path}", shell=True)

    cat_proc = subprocess.Popen(["cat", vcf_path], stdout=subprocess.PIPE)
    annotate_proc = subprocess.Popen([bcftools, "annotate", "-a", bed_path, "-c", "CHROM,-,POS,ID"],
                                     stdin=cat_proc.stdout, stdout=subprocess.PIPE)
    bgzip_proc = subprocess.Popen(["bgzip", "-c"], stdin=annotate_proc.stdout, stdout=subprocess.PIPE)

    with open(out_path, "wb") as f:
        f.write(bgzip_proc.communicate()[0])

    # Generate AFO Geno 38 file
    gen_data = {}
    counter_template = {}

    # with gzip.open(ref_bed, "rt") as IN:
    with subprocess.Popen(["gzip", "-dcf", ref_bed], stdout=subprocess.PIPE).stdout as IN:
        for line in IN:
            line = line.strip()
            if not line:
                continue

            # data = line.split(b"\t")
            data = line.decode().split("\t")
            # gen_data.setdefault(data[0], {})[int(data[2])] = data
            gen_data[(data[0], data[2])] = data
            counter_template[(data[0], data[2])] = 0

    VCFfilename = f"{resultsdir}/temp/{sample}.genotypes.rsIDs.vcf.gz"
    debugCounter = {}

    with subprocess.Popen(
            f'gzip -dcf {VCFfilename} | grep -v "0/0:0:0:0:0,0,0" | grep -v LowQual | cut -f 1,2,4,5,10 | uniq',
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc1, \
            open(f"{resultsdir}/temp/{sample}.afogeno38.txt", "w") as out_file:
        # for line_bytes in proc1.stdout:
        for line_bytes in proc1.stdout:
            line = line_bytes.decode().rstrip()

            if not line or line.startswith("#"):
                continue

            line = line.replace(",<NON_REF>", "")
            data = line.split("\t")

            tmp = data[4].split(":")
            data[4] = tmp[0]

            genotype = ["NA", "NA"]
            extra = "NA"

            pos = (data[0], data[1])

            if data[4] != "./.":
                a1, a2 = data[2], data[3]
                genotype, extra = build_genotype(a1, a2, data[4])
            else:
                # No reads there
                debugCounter[pos] = debugCounter.get(pos, 0) + 1
                continue

            data.extend(genotype)
            data.append(extra)

            # pos = (data[0], data[1])
            if pos in gen_data:
                d2 = gen_data[pos]
                out_data = "\t".join(d2) + "\t" + "\t".join(data) + "\n"
                out_file.write(out_data)
            else:
                debugCounter[pos] = debugCounter.get(pos, 0) + 1

    # Move the original VCF file to a new location
    os.rename(VCFfilename, f"{resultsdir}/{sample}.genotypingVCF.vcf.gz")


if __name__ == '__main__':
    fire.Fire(afogeno_generator_from_VCF)
