#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use File::chdir;


### Processing Needed steps	###

################### parameters

my $dir="/GenomeChronicler/";
my $template_withVEP = "${dir}/templates/reportTemplate_withVEP.tex";
my $template_ohneVEP = "${dir}/templates/reportTemplate_ohneVEP.tex";
my $template = $template_ohneVEP;


#Take BAM from input;
#Take VEP html from input if exists; (if doesn't exist disable VEP processing and choose the appropriate LaTeX template)
#Potentially accept other options like SAMPLE name override


###Sort the long opts command line arguments
#my $conf_file = undef; #"ChroniclerConfig.conf";


#Defining input options and their default values...

#my $conf_file = undef;
my $BAM_file = undef;
my $VEP_file = undef;


my ($helpFlag, $debugFlag);


GetOptions(

	'debug|d'  => \$debugFlag,               # Help/usage
	'help|h'   => \$helpFlag,                # Help/usage

#	'configFile=s' => \$conf_file,
	'bamFile|bam=s' => \$BAM_file,
	'vepFile|vep|html=s' => \$VEP_file,

);



# Get the script filename for use in usage messages
my $scriptName = basename( $0, () );
my $dtag = `date '+%y-%j'`;
chomp($dtag);

#Start timer
my $start_time = time();


##################### check file existence and proceed

if(!defined($BAM_file)) {
	&headerascii();
    &usage();

	print STDERR "\t --- ERROR: No BAM file specified. Please check the usage notes above and try again ---\n";
	exit(500);
}

if(! -e ($BAM_file)) {
	&headerascii();
	print STDERR "\t --- ERROR: The BAM file specified in the command line wasn't found [ $BAM_file ], please check the provided path and try again ---\n";
	exit(404);
}

if(defined($VEP_file) and (!-e ($VEP_file))) {
	&headerascii();
	print STDERR "\t --- ERROR: The VEP file specified in the command line wasn't found [ $VEP_file ], please check the provided path and try again. If you don't want to run the report with VEP, just omit this parameter ---\n";
	exit(606);
}

#wishlist: For V2, make sure the needed software is also in place and callable from the script

#(my $sample= basename($BAM_file)) =~ s/\.[^.]+$//;
#$sample =~ s/\.recal//g;
#system("mkdir -p results/results_${sample}/temp");


(my $sample= basename($BAM_file)) =~ s/\.[^.]+$//;
$sample =~ s/\.recal//g;
$sample =~ s/\.bam\.clean//gi;

system("mkdir -p ${dir}/results/results_${sample}/temp");

#my $LOGFILE1 = "${dir}/results/results_${sample}/${sample}.processingLog.stdout.txt";
my $LOGFILE2 = "${dir}/results/results_${sample}/${sample}.processingLog.stderr.txt";

# Dump parameters
if ( $debugFlag ) {
    print STDERR "SCRIPT = $scriptName\n";
    print STDERR "DEBUG = $debugFlag\n";
    print STDERR "HELP = $helpFlag\n";
    print STDERR "BAM = $BAM_file\n";
    print STDERR "VEP = $VEP_file\n\n";

    print STDERR "SAMPLE = $sample\n";
    print STDERR "OUTDIR = ${dir}/results/results_${sample}/temp\n";
  #  print STDERR "LOGFILE1 = $LOGFILE1\n";
    print STDERR "LOGFILE2 = $LOGFILE2\n\n";

    exit(0);
}




# Print usage and exit if requested
if ( $helpFlag ) {
	&headerascii();
    &usage();
    exit(0);
}



#	--configFile=					MapMi configuration file
#	Note: Specifying a configuration file will override every other option given as input.

sub headerascii {
    print STDERR <<EOP

 #####                                         #####                                                            
#     # ###### #    #  ####  #    # ######    #     # #    # #####   ####  #    # #  ####  #      ###### #####  
#       #      ##   # #    # ##  ## #         #       #    # #    # #    # ##   # # #    # #      #      #    # 
#  #### #####  # #  # #    # # ## # #####     #       ###### #    # #    # # #  # # #      #      #####  #    # 
#     # #      #  # # #    # #    # #         #       #    # #####  #    # #  # # # #      #      #      #####  
#     # #      #   ## #    # #    # #         #     # #    # #   #  #    # #   ## # #    # #      #      #   #  
 #####  ###### #    #  ####  #    # ######     #####  #    # #    #  ####  #    # #  ####  ###### ###### #    #  

\tCopyright(C) 2016-2019 Jose Afonso Guerra-Assuncao et al \@ Personal Genome Project - United Kingdom
\tSee also: https://www.personalgenomes.org.uk

EOP
;
}


sub usage {
    print STDERR <<EOF

\t+++ Welcome to GenomeChronicler - Version 19-291 +++

[USAGE]
$scriptName -h
$scriptName --bamFile QualityRecal_BAMfile.bam [ --vepFile vep_summary_from_WGS_variants.html ]


[PARAMETERS]
	--bamFile=		[REQUIRED] The path to a BAM file that has been preprocessed through 
						markDuplicates and VariantQualityScoreRecalibration. This can be
						obtained by running the first step of the Sarek nextflow pipeline,
						or through other means that do respect the general principles of
						the GATK Variation Calling Best Practices workflow. Note that no
						variation calling is needed to run GenomeChronicler.

	--vepFile=		[OPTIONAL] For the summary tables to appear in the report, a VEP 
						summary HTML file must be provided. This will likely be generated
						if the data is from whole genome sequencing and variants were called
						(e.g. by running all the germline calling steps of the Sarek nextflow
						pipeline or other GATK Best Practices based workflow). If this isn't
						provided, summary tables and plots will automatically be excluded from
						the final report.

	-h, --help		Prints this help page    
  
EOF
;
}


# if(defined($config_file)) {
# 	use Config::General;

# 	#Read config file from disk
# 	my $conf = new Config::General($conf_file);
# 	my %config = $conf->getall;

# 	########## Setup parameters with data from config file
# 	$query_file = $config{"default_query_file"};

# 	#Define rules for handling reads overlapping the loop. Added default values so old configuration files are still valid
# 	$loop_overlap_exclude = $config{"loop_overlap_exclude"} if(defined($config{"loop_overlap_exclude"}));
# 	$loop_overlap_threshold = $config{"loop_overlap_threshold"} if(defined($config{"loop_overlap_threshold"}));

# }




###################################### Main bit of code
&headerascii();

print STDERR "\t +++ INFO: Starting Processing at $start_time\n";
print STDERR "\t +++ INFO: Opening Log File at: $LOGFILE2\n";
system("echo > $LOGFILE2");

my $cleanSample = $sample;
$cleanSample =~ s/\_/\\_/gi;

# 	@echo "${SAMPLE}" | tr -d "\n" > SampleName.txt
open(OUT, ">${dir}/results/results_${sample}/SampleName.txt") or die "Could not open input file: $!\n";
print OUT "$cleanSample";
close OUT;




###################### If we are running VEP, call script or function to process the HTML into the needed tables
if(defined($VEP_file)) {

	$template = $template_withVEP;

	print STDERR "\t +++ INFO: Preprocessing VEP file\n";
	system("perl ${dir}scripts/GenomeChronicler_vepTables_fromVEP.pl $VEP_file results/results_${sample}/");

}


###################### Run script or subroutine to make sure that we don't have 'chr' in the contig names

if(defined($BAM_file)) {

	print STDERR "\t +++ INFO: Preprocessing BAM file\n";	
	
	if((!-e $BAM_file.".clean.BAM") or (!-e $BAM_file.".clean.BAM.bai")) {
		&cleanBAMfile_noCHR($BAM_file);
	}

	$BAM_file = $BAM_file.".clean.BAM";
}


##################### Use the BAM to call the right script to compute currentAncestryPlot

print STDERR "\t +++ INFO: Generating Ancestry\n";

system("perl ${dir}/scripts/GenomeChronicler_ancestry_generator_fromBAM.pl $BAM_file 2>>$LOGFILE2");
system("SAMPLE=$sample ID=$sample DIR=$dir R CMD BATCH ${dir}/scripts/GenomeChronicler_plot_generator_fromAncestry.R");

##################### Use the BAM to call the genotypes on the needed positions for this

print STDERR "\t +++ INFO: Generating Genotypes Files\n";

system("perl ${dir}/scripts/GenomeChronicler_afogeno_generator_fromBAM.pl $BAM_file 2>>$LOGFILE2");
my $AFOgeno_file = "${dir}/results/results_${sample}/temp/${sample}.afogeno38.txt";


##################### Use the generated genotypes file to produce the report tables by linking with the databases

print STDERR "\t +++ INFO: Generating Genome Report Tables\n";

system("perl ${dir}/scripts/GenomeChronicler_genoTables_fromAfoGeno.pl $AFOgeno_file ${dir}/results/results_${sample}/ 2>>$LOGFILE2");


##################### Table filtering for variants that have 0 magnitude and/or are unsupported by external links.


print STDERR "\t +++ INFO: Filtering Report Tables\n";

system("perl ${dir}/scripts/GenomeChronicler_quickFilterFinalReportTables.pl ${dir}/results/results_${sample}/latest.good.reportTable.csv");
system("perl ${dir}/scripts/GenomeChronicler_quickFilterFinalReportTables.pl ${dir}/results/results_${sample}/latest.bad.reportTable.csv");
system("perl ${dir}/scripts/GenomeChronicler_quickFilterFinalReportTables.pl ${dir}/results/results_${sample}/latest.genoset.reportTable.csv");


##################### Call script to summarise found phenotypes as XLS spreadsheet

print STDERR "\t +++ INFO: Combining Excel Tables\n";

system("perl ${dir}/scripts/GenomeChronicler_XLSX_fromTables.pl results/results_${sample}/ ${dir}/results/results_${sample}/${sample}_genotypes_${dtag}.xlsx");


##################### Call LaTeX on the right template to produce the final report

print STDERR "\t +++ INFO: Compiling Genome Report\n";

system("cp $template ${dir}/results/results_${sample}/${sample}_report_${dtag}.tex");
system("cp ${dir}/templates/versionTable.txt ${dir}/results/results_${sample}/");
system("cp ${dir}/templates/GeneStructure.pdf ${dir}/results/results_${sample}/");

my $TEMPLATETEX = "${sample}_report_".$dtag;

&runLatex();

sub runLatex {
    local $CWD = "${dir}/results/results_${sample}";

	for(my $i = 0; $i < 3 ; $i++) {
		system("pdflatex -interaction=nonstopmode ${TEMPLATETEX} .tex 2> /dev/null >/dev/null");
	}
}

##################### Clean up temp files, both in the results folder and in the potentially extraneous tables produced for the report

print STDERR "\t +++ INFO: Cleaning up Temporary and Intermediate Files\n";

system("rm -rf $BAM_file ${BAM_file}.bai");
system("rm -rf ${dir}/results/results_${sample}/temp/");
system("rm -rf ${dir}/results/results_${sample}/latest*.csv");
system("rm -rf ${dir}/results/results_${sample}/versionTable.txt ${dir}/results/results_${sample}/GeneStructure.pdf");
system("rm -rf ${dir}/results/results_${sample}/${TEMPLATETEX}.out ${dir}/results/results_${sample}/texput.log ${dir}/results/results_${sample}/${TEMPLATETEX}.aux ${dir}/results/results_${sample}/${TEMPLATETEX}.log ${dir}/results/results_${sample}/${TEMPLATETEX}.tex");
system("rm -rf ${dir}/GenomeChronicler_plot_generator_fromAncestry.Rout");

sleep (1);

$BAM_file =~ s/\.clean\.BAM//gi;
print STDERR "\n\t +++ DONE: Finished GenomeChronicler for file [ $BAM_file ] in ".( time() - $start_time )." seconds\n";






sub cleanBAMfile_noCHR() {

	my $filename = shift;

	#open(IN, "${dir}/software/samtools view -H $filename |") or die "Could not open input file: $!\n";
	open(IN, "samtools view -H $filename |") or die "Could not open input file: $!\n";
	open(OUT, ">$filename.tempHeader") or die "Could not open output file: $!\n";

	while (my $line = <IN>) {
		chomp($line);
		next if($line eq "");
		
		$line =~ s/SQ\tSN:chr/SQ\tSN:/g;

		print OUT $line."\n";

	}

	close IN;
	close OUT;


	#system("${dir}/software/samtools reheader $filename.tempHeader $filename > $filename.clean.BAM");
	#system("${dir}/software/samtools index $filename.clean.BAM");	
	system("samtools reheader $filename.tempHeader $filename > $filename.clean.BAM");
	system("samtools index $filename.clean.BAM");

	#system("mv $filename.tempBAM $filename");
	system("rm $filename.tempHeader")
}


