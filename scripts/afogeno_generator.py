#!/usr/bin/env python3
import os
import subprocess
import sys
from pathlib import Path

import fire

if __name__ == '__main__':
    sys.path.append(os.path.abspath(os.getcwd()))
from scripts.ancestry_generator import get_gvcf


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


# def uniq(l):
#     return list(set(l))

def prepare_snps_ref_bed(ref_input_path, output_dir, bgzip="bgzip", tabix="tabix"):
    """
    Sort and compress BED file
    - Sort by chromosome and position
    - Compress with bgzip
    - Index with tabix

    Input: .bed
    Output: .bed.gz, .bed.gz.tbi

    Parameters
    ----------
    ref_input_path
    output_dir
    bgzip
    tabix

    Returns
    -------

    """
    ref_input_path = Path(ref_input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ref_path = output_dir / f"{ref_input_path.name}"
    gz_path = output_dir / f"{ref_input_path.name}.gz"
    gz_tbi_path = output_dir / f"{ref_input_path.name}.gz.tbi"

    # Sort and compress BED
    if not gz_path.exists():
        # subprocess.run(["sort", "-k", "1,1", "-k2,2n", str(ref_input_path), ">", str(ref_path)])
        # subprocess.run(["mv", "tempAAA", ref_path])
        # subprocess.run([bgzip, "-c", f"{ref_path}" + " > " + f"{gz_path}"])
        cmd = f"sort -k 1,1 -k2,2n {ref_input_path} > {ref_path}"
        subprocess.run(cmd, shell=True)
        cmd = f"bgzip -c {ref_path} > {gz_path}"
        subprocess.run(cmd, shell=True)

    # Force tabix re-index
    if not gz_tbi_path.exists():
        # tabix index
        # subprocess.run([tabix, "-p", "bed", f"{gz_path}"])
        cmd = f"tabix -p bed {gz_path}"
        subprocess.run(cmd, shell=True)

def generate_afogeno_38(vcf_path, ref_bed, output_dir):
    ## Generate AFO Geno 38 file
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

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # output_vcf_path = output_dir/f"afogeno.genotypes.rsIDs.vcf.gz"
    afogeno_path = output_dir/f"afogeno38.txt"
    debugCounter = {}

    with subprocess.Popen(
            f'gzip -dcf {vcf_path} | grep -v "0/0:0:0:0:0,0,0" | grep -v LowQual | cut -f 1,2,4,5,10 | uniq',
            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as proc1, \
            open(afogeno_path, "w") as out_file:
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


def generate_afogeno_pipeline(bam_or_vcf, output_dir="", numThreads=4,
                              dir_reference='reference', dir_software='software', workdir="",
                              ):
    if 'SINGULARITY_NAME' in os.environ:
        workdir = "/GenomeChronicler/"
    workdir = Path(workdir)

    # resultsdir = workdir
    ref_bed = workdir / f"{dir_reference}/snps.19-114.unique.nochr.bed"
    gatk = workdir / f"{dir_software}/GenomeAnalysisTK.jar"
    ref_hs38 = workdir / f"{dir_reference}/GRCh38_full_analysis_set_plus_decoy_hla_noChr.fa"
    bcftools = "bcftools"
    bgzip = "bgzip"
    tabix = "tabix"

    output_dir = Path(output_dir)
    plink_snps_dir = output_dir / 'prepare_snps_ref_bed'
    prepare_snps_ref_bed(ref_bed, plink_snps_dir, bgzip, tabix)

    ref_snps = plink_snps_dir/'snps.19-114.unique.nochr.bed'
    vcf_dir = output_dir/'vcf_afogeno'
    get_gvcf(bam_or_vcf, vcf_dir, ref_hs38, ref_snps, 'afogeno', gatk, numThreads, bcftools, bgzip)

    ## Generate AFO Geno 38 file
    vcf_afogeno_path = vcf_dir / f'afogeno.rsIDs.gvcf.gz'
    generate_afogeno_38(vcf_afogeno_path, ref_bed, output_dir)


funcs = {
    'prepare_snps_ref_bed': prepare_snps_ref_bed,
    'generate_afogeno_38': generate_afogeno_38,
    'generate_afogeno_pipeline': generate_afogeno_pipeline,
}


class Tools(object):
    def __init__(self):
        super(Tools, self).__init__()


if __name__ == '__main__':
    for k, v in funcs.items():
        setattr(Tools, k, staticmethod(v))
    fire.Fire(Tools)
