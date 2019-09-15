 #####                                         #####                                                            
#     # ###### #    #  ####  #    # ######    #     # #    # #####   ####  #    # #  ####  #      ###### #####  
#       #      ##   # #    # ##  ## #         #       #    # #    # #    # ##   # # #    # #      #      #    # 
#  #### #####  # #  # #    # # ## # #####     #       ###### #    # #    # # #  # # #      #      #####  #    # 
#     # #      #  # # #    # #    # #         #       #    # #####  #    # #  # # # #      #      #      #####  
#     # #      #   ## #    # #    # #         #     # #    # #   #  #    # #   ## # #    # #      #      #   #  
 #####  ###### #    #  ####  #    # ######     #####  #    # #    #  ####  #    # #  ####  ###### ###### #    #  


######## Welcome, thank you for downloading Genome Chronicler, let's do some initial setup ########

### Download reference data

### Uncompress reference data
#tar xvf reference.tar.gz

### Download some needed binaries
wget --http-user=blicp --http-password=bpublicdata http://genomicsdata.cs.ucl.ac.uk/blic_public/afonso/PGPUK/software.tar.gz

### Uncompress software folders
tar xvf software.tar.gz

### Setup symlinks for software (change this to softare.darwin if you are running this on a mac)
ln -sf software.linux software

### Download test data?

### Uncompress test data
#tar xvf testData.tar 

### Run initial test to check for further problems and unmet dependencies
#time perl GenomeChronicler_mainDruid.pl --bamFile=testData/tester_1_GRCh38.bam --vepFile=testData/tester_1_VEP.html