#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from scripts.afogeno_generator import generate_afogeno_pipeline
from scripts.ancestry_generator import generate_ancestry_pipeline
from scripts.genotables_from_afogeno import geno_tables_from_afogeno
from scripts.io import render_latex, export_genotypes_xlsx_from_csvs
from scripts.plot_generator_fromAncestry import plot_generator_from_ancestry
from scripts.utils import clean_bam_file_noCHR
from scripts.vep_tables_from_vep import get_vep_tables_from_vep

def check_path_exist(path, file_type,
                     error_template="--- ERROR: The {file_type} specified in the command line wasn't found [{"
                                    "path}], please check the provided path and try again ---\n",
                     exit_code=404):
    if path and not Path(path).exists():
        print(error_template.format(file_type=file_type, path=path), file=sys.stderr)
        sys.exit(exit_code)


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


def get_parser():
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
    # parser.add_argument('--clean_temporary_files', type=bool, default=True)
    # parser.add_argument('--feature', dest='clean_temporary_files', action='store_true')
    parser.add_argument('--no_clean_temporary_files', default=False, dest='no_clean_temporary_files', action='store_true')
    parser.add_argument('--verbose', dest='verbose', type=bool, default=False)
    return parser

def main():
    import os

    ### Processing Needed steps ###

    ################### parameters
    dir = ""
    if "SINGULARITY_NAME" in os.environ:
        dir = "/GenomeChronicler/"

    resultsdir = os.getcwd()

    # Defining input options and their default values...
    parser = get_parser()
    args = parser.parse_args()

    # conf_file = None
    debugFlag = args.debug
    BAM_file = args.BAM_file
    gVCF_file = args.gVCF_file
    VEP_file = args.VEP_file
    templateParam = args.templateParam
    resultsdir = args.resultsdir
    GATKthreads = args.GATKthreads
    verbose = args.verbose
    no_clean_temporary_files = args.no_clean_temporary_files


    if args.templateParam:
        template = args.templateParam

    scriptName = Path(__file__).name
    dtag = datetime.now().strftime("%y-%j")

    # Start timer
    start_time = datetime.now()

    ##################### check file existence and proceed
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

    output_dir = Path(resultsdir)
    temp_dir = output_dir/"temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # os.system(f"mkdir -p {temp_dir}")
    # check_path_exist(temp_dir, "Results directory", exit_code=102,
    #                  error_template="\t --- ERROR: Results directory [{path}] can't be found and could not be "
    #                                 "created. Please check permissions and try again ---\n")

    LOGFILE2 = f"{output_dir}/{sample}.processingLog.stderr.txt"

    # Main bit of code
    print_header_ascii()
    print(f"\t +++ INFO: Starting Processing at {start_time.strftime('%Y-%m-%d %H:%M:%S')}", file=sys.stderr)
    print(f"\t +++ INFO: Opening Log File at: {LOGFILE2}", file=sys.stderr)

    subprocess.run(f"echo > {LOGFILE2}", shell=True)

    output_vep_dir = None
    if VEP_file:
        print("\t +++ INFO: Preprocessing VEP file")
        output_vep_dir = output_dir/"vep_dir"
        get_vep_tables_from_vep(VEP_file, output_vep_dir)

    if BAM_file:
        print("\t +++ INFO: Preprocessing BAM file")

        output_path = temp_dir/f"{sample}.clean.BAM"
        if not output_path.exists():
            clean_bam_file_noCHR(BAM_file, temp_dir)

        input_path = output_path
    elif gVCF_file:
        print("\t +++ INFO: Preprocessing VCF file")
        input_path = gVCF_file
    else:
        print("\t +++ ERROR: No BAM or VCF file provided. Exiting.")
        sys.exit()

    print("\t +++ INFO: Generating Ancestry")
    output_ancestry_dir = output_dir/"ancestry"
    output_ancestry_path = output_ancestry_dir/'vcf_ancestry'/'ancestry.rsIDs.gvcf.gz'
    if not output_ancestry_path.exists():
        # cmd = f"python3 {dir}scripts/ancestry_generator.py ancestry_generator {BAM_file} {output_dir} {GATKthreads} 2>> {LOGFILE2}"
        # subprocess.run(cmd, shell=True)
        generate_ancestry_pipeline(input_path, output_ancestry_dir, GATKthreads)

    print("\t +++ INFO: Generating Ancestry Plots")
    pca_eigenvec_path = f"{output_ancestry_dir}/pca_ancestry/1kGP_pruned_pca_20.eigenvec"
    output_plot_path = f"{output_dir}/AncestryPlot.pdf"
    if not Path(output_plot_path).exists():
        # cmd = f"PCA_PATH={pca_eigenvec_path} ID={sample} OUTPUT_DIR={output_plot_path} R CMD BATCH {dir}scripts/plot_generator_fromAncestry.R"
        # ret = subprocess.run(cmd, shell=True)
        # print(ret)
        plot_generator_from_ancestry(pca_eigenvec_path, sample, output_dir)
        # exit()

    print("\t +++ INFO: Generating Genotypes Files")
    output_afogeno_dir = output_dir/"afogeno"
    output_afogeno_path = output_ancestry_dir/'vcf_afogeno'/'afogeno.rsIDs.gvcf.gz'
    afogeno_file = output_afogeno_dir/f"afogeno38.txt"
    if not afogeno_file.exists():
        # cmd = f"python3 {dir}scripts/ancestry_generator.py ancestry_generator {BAM_file} {output_dir} {GATKthreads} 2>> {LOGFILE2}"
        # subprocess.run(cmd, shell=True)
        generate_afogeno_pipeline(input_path, output_afogeno_dir, GATKthreads)

    print("\t +++ INFO: Generating Genome Report Tables")
    output_tables_dir = output_dir/"tables"
    geno_tables_from_afogeno(afogeno_file, output_tables_dir)
    for path in output_tables_dir.glob("*.csv"):
        shutil.copy(path, output_dir)

    print("\t +++ INFO: Combining Excel Tables")
    output_xlsx_path = output_tables_dir/f"{sample}.genotypes.{dtag}.xlsx"
    export_genotypes_xlsx_from_csvs(output_tables_dir, output_xlsx_path)
    shutil.copy(output_xlsx_path, output_dir) # copy to main results dir

    # render latex template
    print("\t +++ INFO: Compiling Genome Report")
    template_withVEP = f"{dir}templates/reportTemplate_withVEP.tex"
    template_ohneVEP = f"{dir}templates/reportTemplate_ohneVEP.tex"
    template = template_ohneVEP
    if VEP_file:
        template = template_withVEP if not templateParam else templateParam
    run_latex_dir = output_dir/"latex"
    render_latex(output_tables_dir, sample, run_latex_dir,
                 plot_path=output_plot_path, vep_dir=output_vep_dir,
                 template_dir="templates", latex_template=template,
                 )
    latex_pdf_path = run_latex_dir/f"{sample}_report.pdf"
    shutil.copy(latex_pdf_path, output_dir) # copy to main results dir


    # if no_clean_temporary_files is False:
    #     print("\t +++ INFO: Cleaning up Temporary and Intermediate Files")
    #     if BAM_file is not None:
    #         os.remove(BAM_file)
    #         os.remove(f"{BAM_file}.bai")
    #
    #     shutil.rmtree(f"{temp_dir}")
    #     for file in Path(f"{output_dir}/").glob("latest*.csv"):
    #         os.remove(file)
    #
    #     subprocess.run(f"rm -rf {output_dir}/versionTable.txt", shell=True)
    #     subprocess.run(f"rm -rf {output_dir}/GeneStructure.pdf", shell=True)
    #     subprocess.run(f"rm -rf {output_dir}/{TEMPLATETEX}.out", shell=True)
    #     subprocess.run(f"rm -rf {output_dir}/texput.log", shell=True)
    #     subprocess.run(f"rm -rf {output_dir}/{TEMPLATETEX}.aux", shell=True)
    #     subprocess.run(f"rm -rf {output_dir}/{TEMPLATETEX}.log", shell=True)
    #     subprocess.run(f"rm -rf {output_dir}/{TEMPLATETEX}.tex", shell=True)
    #     subprocess.run(f"rm -rf {dir}GenomeChronicler_plot_generator_fromAncestry.Rout", shell=True)

    # # time.sleep(1)
    # if BAM_file is not None:
    #     BAM_file = BAM_file.replace(".clean.BAM", "")
    #     print(
    #         f"\n\t +++ DONE: Finished GenomeChronicler for file [ {BAM_file} ] in {datetime.now() - start_time} seconds")


if __name__ == '__main__':
    main()

