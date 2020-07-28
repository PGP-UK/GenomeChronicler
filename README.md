	 #####                                         #####                                                            
	#     # ###### #    #  ####  #    # ######    #     # #    # #####   ####  #    # #  ####  #      ###### #####  
	#       #      ##   # #    # ##  ## #         #       #    # #    # #    # ##   # # #    # #      #      #    # 
	#  #### #####  # #  # #    # # ## # #####     #       ###### #    # #    # # #  # # #      #      #####  #    # 
	#     # #      #  # # #    # #    # #         #       #    # #####  #    # #  # # # #      #      #      #####  
	#     # #      #   ## #    # #    # #         #     # #    # #   #  #    # #   ## # #    # #      #      #   #  
	 #####  ###### #    #  ####  #    # ######     #####  #    # #    #  ####  #    # #  ####  ###### ###### #    #  




[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/3664)


# Welcome

This is the repository for Genome Chronicler, the Personal Genome Project United Kingdom (PGP-UK) genomic report generation scripts. 

# Getting Started

After cloning this repository, run the SetupMeFirst.sh script in your local system to retrieve the extra files needed to run the pipeline (around 10GB, so too big for git). 

# Input files

The main script (GenomeChronicler_mainDruid.pl) needs a BAM file as input, and optionally can also use a VEP generated summary html file, if variants have already been called on the data and summaries are to be produced.

# Dependencies

To handle the myriad dependencies present in this pipeline, it is avaliable through Singularity Hub as a singularity container (see badge at top of the page). 


# Easy Start using Singularity

If you don't already have singularity on your system, or want to know more about it, head to their userguide at: https://sylabs.io/guides/3.1/user-guide/
While Singularity is not needed to run GenomeChronicler, it does make setup much easier. 

For a manual installation without Singularity, please follow the steps in the %post section of the Singularity file in this repository, to install all the dependencies.



Downloading pre-packaged GenomeChronicler from SingularityHub
````
singularity pull shub://PGP-UK/GenomeChronicler
````


Getting some test data (NA12878 from ENA, pre-mapped to GRCh38, and the respective reference)
````wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram #A bit too large for a starting test
singularity exec GenomeChronicler_latest.sif wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram

singularity exec GenomeChronicler_latest.sif wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

````

Converting data to BAM format 
````
singularity exec GenomeChronicler_latest.sif samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa -b -o NA12878wxs.bam NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram
````


Running GenomeChronicler on the data
````
singularity run GenomeChronicler_latest.sif --bamFile=NA12878wxs.bam 
````

# Command Line Options



|       Option      | Requirement | Description                                                                                                                                                                                                                                                                                                                                                                                                                      |
|:-----------------:|:------------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --bamFile         | REQUIRED     | The path to a BAM file that has been preprocessed through markDuplicates and VariantQualityScoreRecalibration. This can be obtained by running the first step of the Sarek nextflow pipeline, or through other means that do respect the general principles of the GATK Variation Calling Best Practices workflow. Note that no variation calling is needed to run GenomeChronicler.                                             |
| --vepFile         | OPTIONAL     | For the summary tables to appear in the report, a VEP summary HTML file must be provided. This will likely be generated if the data is from whole genome sequencing and variants were called (e.g. by running all the germline calling steps of the Sarek nextflow pipeline or other GATK Best Practices based workflow). If this isn't provided, summary tables and plots will automatically be excluded from the final report. |
| --resultsDir      | OPTIONAL     | For setting the absolute path of the results folder to be produced when running GenomeChronicler.                                                                                                                                                                                                                                                                                                                                |
| --customTemplate  | OPTIONAL     | For customising the output report, set this variable to the path of a custom LaTeX file to act as a template for the report. The default templates bundled with this software can also be found in the project github page.                                                                                                                                                                                                      |
| --GATKthreads     | OPTIONAL     | Number of threads to use for the GATK genotyping steps of this processing pipeline.                                                                                                                                                                                                                                                                                                                                              |
