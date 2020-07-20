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

Downloading pre-packaged GenomeChronicler from SingularityHub
````
singularity pull shub://PGP-UK/GenomeChronicler
````


Getting some test data (NA12878 from ENA, pre-mapped to GRCh38)
````wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram #ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR194/ERR194160/NA12891_S1.bam
singularity exec GenomeChronicler_latest.sif wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram

````

Running GenomeChronicler on the data
````
singularity exec GenomeChronicler_latest.sif samtools -b -o NA12878wxs.bam NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram
singularity run --app gc GenomeChronicler_latest.sif --bamFile=NA12878wxs.bam 
````

# 
