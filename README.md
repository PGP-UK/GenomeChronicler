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

