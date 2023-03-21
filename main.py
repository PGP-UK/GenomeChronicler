#!/usr/bin/env python3

import argparse
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
import shutil
import time
import os


def print_header_ascii():
    print("""
 #####                                         #####                                                            
#     # ###### #    #  ####  #    # ######    #     # #    # #####   ####  #    # #  ####  #      ###### #####  
#       #      ##   # #    # ##  ## #         #       #    # #    # #    # ##   # # #    # #      #      #    # 
#  #### #####  # #  # #    # # ## # #####     #       ###### #    # #    # # #  # # #      #      #####  #    # 
#     # #      #  # # #    # #    # #         #       #    # #####  #    # #  # # # #      #      #      #####  
#     # #      #   ## #    # #    # #         #     # #    # #   #  #    # #   ## # #    # #      #      #   #  
 #####  ###### #    #  ####  #    # ######     #####  #    # #    #  ####  #    # #  ####  ###### ###### #    #  

\tCopyright(C) 2016-2022 Jose Afonso Guerra-Assuncao et al @ Personal Genome Project - United Kingdom
\tSee also: https://www.personalgenomes.org.uk
""")


def cleanBAMfile_noCHR(filename):
    # Create a temporary header file
    with open(f"{filename}.tempHeader", "w") as temp_header_file:
        # Run samtools view command and get its output
        cmd = f"samtools view -H {filename}"
        samtools_view_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, shell=True)

        for line in samtools_view_process.stdout:
            line = line.strip()
            if line == "":
                continue

            # Replace 'SQ\tSN:chr' with 'SQ\tSN:'
            line = line.replace("SQ\tSN:chr", "SQ\tSN:")

            # Write to the temporary header file
            temp_header_file.write(line + "\n")

        samtools_view_process.communicate()

    # Run samtools reheader command
    cmd = f"samtools reheader {filename}.tempHeader {filename}"
    subprocess.run(cmd, stdout=subprocess.PIPE, text=True, shell=True)

    # Create the clean BAM file
    cmd = f"samtools index {filename}.clean.BAM"
    subprocess.run(cmd, shell=True)

    # Remove the temporary header file
    cmd = f"rm {filename}.tempHeader"
    subprocess.run(cmd, shell=True)


def main_druid():
    import os

    ### Processing Needed steps ###

    ################### parameters
    dir = ""
    if "SINGULARITY_NAME" in os.environ:
        dir = "/GenomeChronicler/"

    resultsdir = os.getcwd()
    template_withVEP = f"{dir}templates/reportTemplate_withVEP.tex"
    template_ohneVEP = f"{dir}templates/reportTemplate_ohneVEP.tex"
    template = template_ohneVEP

    # Defining input options and their default values...

    parser = argparse.ArgumentParser(description="GenomeChronicler")
    parser.add_argument('-d', '--debug', action='store_true', help='Debug flag')
    # parser.add_argument('-h', '--help', action='help', help='Help/usage')
    parser.add_argument('--bamFile', '--bam', dest='BAM_file', type=str, help="""
                        [REQUIRED] The path to a BAM file that has been preprocessed through 
                        markDuplicates and VariantQualityScoreRecalibration. This can be
                        obtained by running the first step of the Sarek nextflow pipeline,
                        or through other means that do respect the general principles of
                        the GATK Variation Calling Best Practices workflow. Note that no
                        variation calling is needed to run GenomeChronicler.	
                        """)
    parser.add_argument('--vcfFile', '--vcf', '--gvcf', dest='gVCF_file', type=str, help="""
                        [REQUIRED] The path to a gVCF file produced by GATK or the GRCh38
                        reference genome. This can be obtained by running all steps from 
                        the Sarek nextflow pipeline, or through other means that do 
                        respect the general principles of the GATK Variation Calling Best 
                        Practices workflow. This avoids the need to run the GATK with GC.""")
    parser.add_argument('--vepFile', '--vep', '--html', dest='VEP_file', type=str, help="""
                        [OPTIONAL] For the summary tables to appear in the report, a VEP 
                        summary HTML file must be provided. This will likely be generated
                        if the data is from whole genome sequencing and variants were called
                        (e.g. by running all the germline calling steps of the Sarek nextflow
                        pipeline or other GATK Best Practices based workflow). If this isn't
                        provided, summary tables and plots will automatically be excluded from
                        the final report.""")
    parser.add_argument('--customTemplate', '--template', '--latex', dest='templateParam', type=str,
                        help="""
                        [OPTIONAL] For customising the output report, set this variable to the
                        path of a custom LaTeX file to act as a template for the report. The 
                        default templates bundled with this software can also be found in the 
                        project github page.
                        """)
    parser.add_argument('--resultsDir', '--outputDir', dest='resultsdir', type=str, help="""
                        [OPTIONAL] For setting the absolute path of the results folder to be 
                        produced when running GenomeChronicler.""")
    parser.add_argument('--threads', '--GATKthreads', '-t', dest='GATKthreads', type=int, default=1, help="""
                        [OPTIONAL] Number of threads to use for the GATK genotyping steps of this
                        processing pipeline.""")
    parser.add_argument('--clean_temporary_files', type=bool, default=True)
    # parser.add_argument('--feature', dest='clean_temporary_files', action='store_true')
    parser.add_argument('--no_clean_temporary_files', dest='clean_temporary_files', action='store_false')

    args = parser.parse_args()
    # conf_file = None
    debugFlag = args.debug
    BAM_file = args.BAM_file
    gVCF_file = args.gVCF_file
    VEP_file = args.VEP_file
    templateParam = args.templateParam
    resultsdir = args.resultsdir
    GATKthreads = args.GATKthreads
    clean_temporary_files = args.clean_temporary_files

    if args.templateParam:
        template = args.templateParam

    scriptName = Path(__file__).name
    dtag = datetime.now().strftime("%y-%j")

    # Start timer
    start_time = datetime.now()

    ##################### check file existence and proceed

    def check_path_exist(path, file_type,
                         error_template="--- ERROR: The {file_type} specified in the command line wasn't found [{"
                                        "path}], please check the provided path and try again ---\n",
                         exit_code=404):
        if path and not Path(path).exists():
            print(error_template.format(file_type=file_type, path=path), file=sys.stderr)
            sys.exit(exit_code)

    if BAM_file:
        check_path_exist(BAM_file, "BAM file", exit_code=404)
    elif gVCF_file:
        check_path_exist(gVCF_file, "gVCF file", exit_code=404)
    elif not BAM_file and not gVCF_file:
        print_header_ascii()
        parser.print_help()
        print("\t --- ERROR: No BAM or gVCF file specified. Please check the usage notes above and try again ---\n",
              file=sys.stderr)
        sys.exit(500)

    check_path_exist(VEP_file, "VEP file", exit_code=606,
                     error_template="\t --- ERROR: The {fill_type} specified in the command line wasn't found [{"
                                    "path}], please check the provided path and try again. If you don't want to "
                                    "run the report with VEP, just omit this parameter ---\n")
    check_path_exist(templateParam, "LaTeX Template file", exit_code=501)

    sample = None
    import os
    if BAM_file:
        sample = Path(BAM_file).name.split('.')[0]
        # sample = os.path.splitext(os.path.basename(BAM_file))[0]
        # sample = re.sub(r'\.recal', '', sample)
        # sample = re.sub(r'\.bam\.clean', '', sample, flags=re.IGNORECASE)
    elif gVCF_file:
        sample = Path(gVCF_file).name.split('.')[0]
        # sample = os.path.splitext(os.path.basename(gVCF_file))[0]
        # sample = re.sub(r'\.g\.', '.', sample, flags=re.IGNORECASE)
        # sample = re.sub(r'\.g\.vcf', '', sample, flags=re.IGNORECASE)
        # sample = re.sub(r'\.vcf', '', sample, flags=re.IGNORECASE)
    else:
        print_header_ascii()
        print("\t --- MAJOR ERROR: No BAM or gVCF file found. Please check the usage notes above and try again ---\n",
              file=sys.stderr)
        sys.exit(555)

    output_sample_dir = f"{resultsdir}/results/results_{sample}"
    temp_dir = f"{output_sample_dir}/temp"

    os.system(f"mkdir -p {temp_dir}")
    check_path_exist(temp_dir, "Results directory", exit_code=102,
                     error_template="\t --- ERROR: Results directory [{path}] can't be found and could not be "
                                    "created. Please check permissions and try again ---\n")

    LOGFILE2 = f"{output_sample_dir}/{sample}.processingLog.stderr.txt"

    # Dump parameters
    if debugFlag:
        print(f"SCRIPT = {scriptName}", file=sys.stderr)
        print(f"DEBUG = {debugFlag}", file=sys.stderr)
        print(f"BAM = {BAM_file}", file=sys.stderr)
        print(f"VCF = {gVCF_file}", file=sys.stderr)
        print(f"VEP = {VEP_file}\n", file=sys.stderr)
        print(f"GATKthreads = {GATKthreads}\n\n", file=sys.stderr)

        print(f"TEMPLATE = {template}", file=sys.stderr)
        print(f"SAMPLE = {sample}", file=sys.stderr)
        print(f"DIR = {dir}", file=sys.stderr)
        print(f"TMPDIR = {temp_dir}\n", file=sys.stderr)
        #  print(f"LOGFILE1 = {LOGFILE1}\n", file=sys.stderr)
        print(f"LOGFILE2 = {LOGFILE2}\n\n", file=sys.stderr)

        sys.exit(0)

    # Main bit of code
    print_header_ascii()

    print(f"\t +++ INFO: Starting Processing at {start_time.strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stderr)
    print(f"\t +++ INFO: Opening Log File at: {LOGFILE2}", file=sys.stderr)

    subprocess.run(f"echo > {LOGFILE2}", shell=True)

    clean_sample = sample.replace("_", "\\_")

    with open(f"{output_sample_dir}/SampleName.txt", "w") as sample_name_file:
        sample_name_file.write(clean_sample)

    if VEP_file:
        template = template_withVEP if not templateParam else templateParam

        print("\t +++ INFO: Preprocessing VEP file")
        cmd = f"python3 {dir}scripts/GenomeChronicler_vepTables_fromVEP.py {VEP_file} {output_sample_dir}/"
        subprocess.run(cmd, shell=True)

    if BAM_file:
        print("\t +++ INFO: Preprocessing BAM file")

        if not os.path.exists(f"{BAM_file}.clean.BAM") or not os.path.exists(f"{BAM_file}.clean.BAM.bai"):
            cleanBAMfile_noCHR(BAM_file)

        BAM_file = f"{BAM_file}.clean.BAM"

        print("\t +++ INFO: Generating Ancestry")
        cmd = f"python3 {dir}scripts/GenomeChronicler_ancestry_generator_fromBAM.py {BAM_file} {output_sample_dir}/ {GATKthreads} 2>> {LOGFILE2}"
        subprocess.run(cmd, shell=True)
        cmd = f"SAMPLE={sample} ID={sample} DIR={resultsdir} R CMD BATCH {dir}scripts/GenomeChronicler_plot_generator_fromAncestry.R"
        subprocess.run(cmd, shell=True)

        print("\t +++ INFO: Generating Genotypes Files")
        cmd = f"python3 {dir}scripts/GenomeChronicler_afogeno_generator_fromBAM.py {BAM_file} {output_sample_dir}/ {GATKthreads} 2>> {LOGFILE2}"
        subprocess.run(cmd, shell=True)
        AFOgeno_file = f"{temp_dir}/{sample}.afogeno38.txt"

    elif gVCF_file:
        print("\t +++ INFO: Preprocessing VCF file")

        print("\t +++ INFO: Generating Ancestry")
        if not Path(f"{output_sample_dir}/{sample}.genotypingVCF.vcf.gz").exists():
            cmd = f"python3 {dir}scripts/GenomeChronicler_ancestry_generator_fromVCF.py {gVCF_file} {sample} {output_sample_dir} 2>> {LOGFILE2}"
            subprocess.run(cmd, shell=True)
        # if not Path(f"{output_sample_dir}/{sample}.genotypingVCF.vcf.gz").exists():
        cmd = f"python3 {dir}scripts/GenomeChronicler_plot_generator_fromAncestry.py {output_sample_dir}/temp/{sample}_1kGP_pruned_pca_20.eigenvec {sample} {output_sample_dir}/"
        subprocess.run(cmd, shell=True)

        print("\t +++ INFO: Generating Genotypes Files")
        cmd = f"python3 {dir}scripts/GenomeChronicler_afogeno_generator_fromVCF.py {gVCF_file} {sample} {output_sample_dir} 2>> {LOGFILE2}"
        subprocess.run(cmd, shell=True)
        AFOgeno_file = f"{temp_dir}/{sample}.afogeno38.txt"

    else:
        print("\t +++ ERROR: No BAM or VCF file provided. Exiting.")
        sys.exit()

    print("\t +++ INFO: Generating Genome Report Tables")
    if not Path(f"{output_sample_dir}/latest.good.reportTable.csv").exists():
        cmd = f"python3 {dir}scripts/GenomeChronicler_genoTables_fromAfoGeno.py {AFOgeno_file} {output_sample_dir}/ 2>> {LOGFILE2}"
        subprocess.run(cmd, shell=True)

    print("\t +++ INFO: Filtering Report Tables")
    cmd = f"python3 {dir}scripts/GenomeChronicler_quickFilterFinalReportTables.py {output_sample_dir}/latest.good.reportTable.csv"
    subprocess.run(cmd, shell=True)
    cmd = f"python3 {dir}scripts/GenomeChronicler_quickFilterFinalReportTables.py {output_sample_dir}/latest.bad.reportTable.csv"
    subprocess.run(cmd, shell=True)

    print("\t +++ INFO: Combining Excel Tables")
    cmd = f"python3 {dir}scripts/GenomeChronicler_XLSX_fromTables.py {output_sample_dir}/ {output_sample_dir}/{sample}genotypes{dtag}.xlsx"
    subprocess.run(cmd, shell=True)

    print("\t +++ INFO: Compiling Genome Report")
    TEMPLATETEX = f"{sample}_report_{dtag}"
    shutil.copy(template, f"{output_sample_dir}/{TEMPLATETEX}.tex")
    shutil.copy(f"{dir}templates/versionTable.txt", f"{output_sample_dir}/")
    shutil.copy(f"{dir}templates/GeneStructure.pdf", f"{output_sample_dir}/")

    def run_latex():
        cwd = os.getcwd()
        os.chdir(f"{resultsdir}/results/results_{sample}")
        for _ in range(3):
            cmd = f"pdflatex -interaction=nonstopmode {TEMPLATETEX}.tex 2> /dev/null > /dev/null"
            subprocess.run(cmd, shell=True)
        os.chdir(cwd)

    run_latex()


    if clean_temporary_files is True:
        print("\t +++ INFO: Cleaning up Temporary and Intermediate Files")
        if BAM_file is not None:
            os.remove(BAM_file)
            os.remove(f"{BAM_file}.bai")

        shutil.rmtree(f"{temp_dir}")
        for file in Path(f"{output_sample_dir}/").glob("latest*.csv"):
            os.remove(file)

        subprocess.run(f"rm -rf {output_sample_dir}/versionTable.txt", shell=True)
        subprocess.run(f"rm -rf {output_sample_dir}/GeneStructure.pdf", shell=True)
        subprocess.run(f"rm -rf {output_sample_dir}/{TEMPLATETEX}.out", shell=True)
        subprocess.run(f"rm -rf {output_sample_dir}/texput.log", shell=True)
        subprocess.run(f"rm -rf {output_sample_dir}/{TEMPLATETEX}.aux", shell=True)
        subprocess.run(f"rm -rf {output_sample_dir}/{TEMPLATETEX}.log", shell=True)
        subprocess.run(f"rm -rf {output_sample_dir}/{TEMPLATETEX}.tex", shell=True)
        subprocess.run(f"rm -rf {dir}GenomeChronicler_plot_generator_fromAncestry.Rout", shell=True)

    time.sleep(1)
    if BAM_file is not None:
        BAM_file = BAM_file.replace(".clean.BAM", "")
        print(
            f"\n\t +++ DONE: Finished GenomeChronicler for file [ {BAM_file} ] in {datetime.now() - start_time} seconds")


if __name__ == '__main__':
    main_druid()
