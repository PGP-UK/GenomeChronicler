#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use File::chdir;


### Processing Needed steps	###

################### parameters

my $dir="";
if(defined $ENV{'SINGULARITY_NAME'}){
	$dir="/GenomeChronicler/";
}


my $resultsdir=`pwd`;
chomp($resultsdir);
my $template_withVEP = "${dir}templates/reportTemplate_withVEP.tex";
my $template_ohneVEP = "${dir}templates/reportTemplate_ohneVEP.tex";
my $template = $template_ohneVEP;
my $GATKthreads = 1;

#Take BAM from input;
#Take VEP html from input if exists; (if doesn't exist disable VEP processing and choose the appropriate LaTeX template)
#Potentially accept other options like SAMPLE name override


###Sort the long opts command line arguments
#my $conf_file = undef; #"ChroniclerConfig.conf";


#Defining input options and their default values...

#my $conf_file = undef;
my $BAM_file = undef;
my $gVCF_file = undef;
my $VEP_file = undef;


my ($helpFlag, $debugFlag, $templateParam);


GetOptions(

	'debug|d'  => \$debugFlag,               # Help/usage
	'help|h'   => \$helpFlag,                # Help/usage

#	'configFile=s' => \$conf_file,
	'bamFile|bam=s' => \$BAM_file,
	'vcfFile|vcf|gvcf=s' => \$gVCF_file,
	'vepFile|vep|html=s' => \$VEP_file,
	'customTemplate|template|latex=s' => \$templateParam,
	'resultsDir|outputDir=s' => \$resultsdir,
	'threads|GATKthreads|t=i' => \$GATKthreads,
);

if(defined($templateParam)) {
	$template = $templateParam;
}

# Get the script filename for use in usage messages
my $scriptName = basename( $0, () );
my $dtag = `date '+%y-%j'`;
chomp($dtag);

#Start timer
my $start_time = time();


##################### check file existence and proceed

if(!defined($BAM_file) and !defined($gVCF_file)) {
	&headerascii();
    &usage();

	print STDERR "\t --- ERROR: No BAM or gVCF file specified. Please check the usage notes above and try again ---\n";
	exit(500);
}

if(defined($BAM_file) and (!-e ($BAM_file))) {
	&headerascii();
	print STDERR "\t --- ERROR: The BAM file specified in the command line wasn't found [ $BAM_file ], please check the provided path and try again ---\n";
	exit(404);
}

if(defined($gVCF_file) and (!-e ($gVCF_file))) {
	&headerascii();
	print STDERR "\t --- ERROR: The gVCF file specified in the command line wasn't found [ $gVCF_file ], please check the provided path and try again ---\n";
	exit(404);
}

if(defined($VEP_file) and (!-e ($VEP_file))) {
	&headerascii();
	print STDERR "\t --- ERROR: The VEP file specified in the command line wasn't found [ $VEP_file ], please check the provided path and try again. If you don't want to run the report with VEP, just omit this parameter ---\n";
	exit(606);
}

if(defined($templateParam) and (!-e ($templateParam))) {
	&headerascii();
	print STDERR "\t --- ERROR: The LaTeX Template file specified in the command line wasn't found [ $templateParam ], please check the provided path and try again. ---\n";
	exit(501);
}

my $sample = undef;
if(defined($BAM_file)) {
	($sample= basename($BAM_file)) =~ s/\.[^.]+$//;
	$sample =~ s/\.recal//g;
	$sample =~ s/\.bam\.clean//gi;
}
elsif(defined($gVCF_file)) {
	($sample= basename($gVCF_file)) =~ s/\.[^.]+$//;
	$sample =~ s/\.g.vcf//gi;
	$sample =~ s/\.vcf//gi;
}
else {
	&headerascii();
	print STDERR "\t --- MAJOR ERROR: No BAM or gVCF file found. Please check the usage notes above and try again ---\n";
	exit(555);
}

system("mkdir -p ${resultsdir}/results/results_${sample}/temp");

if(! -e ($resultsdir)) {
	&headerascii();
	print STDERR "\t --- ERROR: Results directory [$resultsdir] can't be found and could not be created. Please check permissions and try again ---\n";
	exit(102);
}

my $LOGFILE2 = "${resultsdir}/results/results_${sample}/${sample}.processingLog.stderr.txt";

# Dump parameters
if ( $debugFlag ) {
    print STDERR "SCRIPT = $scriptName\n";
    print STDERR "DEBUG = $debugFlag\n";
    print STDERR "HELP = $helpFlag\n";
    print STDERR "BAM = $BAM_file\n";
	print STDERR "VCF = $gVCF_file\n";
    print STDERR "VEP = $VEP_file\n\n";
    print STDERR "GATKthreads = $GATKthreads\n\n";

    print STDERR "TEMPLATE = $template\n";
    print STDERR "SAMPLE = $sample\n";
    print STDERR "DIR = ${dir}\n";
    print STDERR "TMPDIR = ${resultsdir}/results/results_${sample}/temp\n";
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



sub headerascii {
    print STDERR <<EOP

 #####                                         #####                                                            
#     # ###### #    #  ####  #    # ######    #     # #    # #####   ####  #    # #  ####  #      ###### #####  
#       #      ##   # #    # ##  ## #         #       #    # #    # #    # ##   # # #    # #      #      #    # 
#  #### #####  # #  # #    # # ## # #####     #       ###### #    # #    # # #  # # #      #      #####  #    # 
#     # #      #  # # #    # #    # #         #       #    # #####  #    # #  # # # #      #      #      #####  
#     # #      #   ## #    # #    # #         #     # #    # #   #  #    # #   ## # #    # #      #      #   #  
 #####  ###### #    #  ####  #    # ######     #####  #    # #    #  ####  #    # #  ####  ###### ###### #    #  

\tCopyright(C) 2016-2022 Jose Afonso Guerra-Assuncao et al \@ Personal Genome Project - United Kingdom
\tSee also: https://www.personalgenomes.org.uk

EOP
;
}


sub usage {
    print STDERR <<EOF

\t+++ Welcome to GenomeChronicler - Version 22-146 +++

[USAGE]
$scriptName -h
$scriptName --bamFile QualityRecal_BAMfile.bam [ --vepFile vep_summary_from_WGS_variants.html ]
$scriptName --vcfFile QualityRecal_gVCF.vcf [ --vepFile vep_summary_from_WGS_variants.html ]


[PARAMETERS]
	--bamFile=		[REQUIRED] The path to a BAM file that has been preprocessed through 
						markDuplicates and VariantQualityScoreRecalibration. This can be
						obtained by running the first step of the Sarek nextflow pipeline,
						or through other means that do respect the general principles of
						the GATK Variation Calling Best Practices workflow. Note that no
						variation calling is needed to run GenomeChronicler.	
	--vcfFile=		[REQUIRED] The path to a gVCF file produced by GATK or the GRCh38
						reference genome. This can be obtained by running all steps from 
						the Sarek nextflow pipeline, or through other means that do 
						respect the general principles of the GATK Variation Calling Best 
						Practices workflow. This avoids the need to run the GATK with GC.

	--vepFile=		[OPTIONAL] For the summary tables to appear in the report, a VEP 
						summary HTML file must be provided. This will likely be generated
						if the data is from whole genome sequencing and variants were called
						(e.g. by running all the germline calling steps of the Sarek nextflow
						pipeline or other GATK Best Practices based workflow). If this isn't
						provided, summary tables and plots will automatically be excluded from
						the final report.

	--resultsDir=		[OPTIONAL] For setting the absolute path of the results folder to be 
						produced when running GenomeChronicler.

	--customTemplate=	[OPTIONAL] For customising the output report, set this variable to the
						path of a custom LaTeX file to act as a template for the report. The 
						default templates bundled with this software can also be found in the 
						project github page.

	--GATKthreads=	[OPTIONAL] Number of threads to use for the GATK genotyping steps of this
						processing pipeline.

	-d, --debug		Prints the relevant variables that are in place after parameter parsing.

	-h, --help		Prints this help page.
  

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
open(OUT, ">${resultsdir}/results/results_${sample}/SampleName.txt") or die "Could not open input file: $!\n";
print OUT "$cleanSample";
close OUT;




###################### If we are running VEP, call script or function to process the HTML into the needed tables
if(defined($VEP_file)) {

	$template = $template_withVEP if(!defined($templateParam));

	print STDERR "\t +++ INFO: Preprocessing VEP file\n";
	system("perl ${dir}scripts/GenomeChronicler_vepTables_fromVEP.pl $VEP_file ${resultsdir}/results/results_${sample}/");

}


###################### Run script or subroutine to make sure that we don't have 'chr' in the contig names
my $AFOgeno_file = undef;

if(defined($BAM_file)) {

	print STDERR "\t +++ INFO: Preprocessing BAM file\n";	
	
	if((!-e $BAM_file.".clean.BAM") or (!-e $BAM_file.".clean.BAM.bai")) {
		&cleanBAMfile_noCHR($BAM_file);
	}

	$BAM_file = $BAM_file.".clean.BAM";



	##################### Use the BAM to call the right script to compute currentAncestryPlot

	print STDERR "\t +++ INFO: Generating Ancestry\n";

#	system("perl ${dir}scripts/GenomeChronicler_ancestry_generator_fromBAM.pl $BAM_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	system("python3 ${dir}scripts/GenomeChronicler_ancestry_generator_fromBAM.py $BAM_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	system("SAMPLE=$sample ID=$sample DIR=$resultsdir R CMD BATCH ${dir}scripts/GenomeChronicler_plot_generator_fromAncestry.R");

	##################### Use the BAM to call the genotypes on the needed positions for this

	print STDERR "\t +++ INFO: Generating Genotypes Files\n";

#	system("perl ${dir}scripts/GenomeChronicler_afogeno_generator_fromBAM.pl $BAM_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	system("python3 ${dir}scripts/GenomeChronicler_afogeno_generator_fromBAM.py $BAM_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	$AFOgeno_file = "${resultsdir}/results/results_${sample}/temp/${sample}.afogeno38.txt";

}
elsif(defined($gVCF_file)) {

	print STDERR "\t +++ INFO: Preprocessing VCF file\n";
	
	# if((!-e $gVCF_file.".clean.vcf") or (!-e $gVCF_file.".clean.vcf.idx")) {
	# 	&cleanVCFfile($gVCF_file);
	# }

	# $gVCF_file = $gVCF_file.".clean.vcf";

	print STDERR "\t +++ INFO: Generating Ancestry\n";

	# system("perl ${dir}scripts/GenomeChronicler_ancestry_generator_fromVCF.pl $gVCF_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	system("python3 ${dir}scripts/GenomeChronicler_ancestry_generator_fromVCF.py $gVCF_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	system("SAMPLE=$sample ID=$sample DIR=$resultsdir R CMD BATCH ${dir}scripts/GenomeChronicler_plot_generator_fromAncestry.R");

	print STDERR "\t +++ INFO: Generating Genotypes Files\n";

	# system("perl ${dir}scripts/GenomeChronicler_afogeno_generator_fromVCF.pl $gVCF_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	system("python3 ${dir}scripts/GenomeChronicler_afogeno_generator_fromVCF.py $gVCF_file $resultsdir $GATKthreads 2>>$LOGFILE2");
	$AFOgeno_file = "${resultsdir}/results/results_${sample}/temp/${sample}.afogeno38.txt";

}
else {
	print STDERR "\t +++ ERROR: No BAM or VCF file provided. Exiting.\n";
	exit;
}


##################### Use the generated genotypes file to produce the report tables by linking with the databases

print STDERR "\t +++ INFO: Generating Genome Report Tables\n";

# system("perl ${dir}scripts/GenomeChronicler_genoTables_fromAfoGeno.pl $AFOgeno_file ${resultsdir}/results/results_${sample}/ 2>>$LOGFILE2");
system("python3 ${dir}scripts/GenomeChronicler_genoTables_fromAfoGeno.py $AFOgeno_file ${resultsdir}/results/results_${sample}/ 2>>$LOGFILE2");


##################### Table filtering for variants that have 0 magnitude and/or are unsupported by external links.


print STDERR "\t +++ INFO: Filtering Report Tables\n";

system("python3 ${dir}scripts/GenomeChronicler_quickFilterFinalReportTables.py ${resultsdir}/results/results_${sample}/latest.good.reportTable.csv");
system("python3 ${dir}scripts/GenomeChronicler_quickFilterFinalReportTables.py ${resultsdir}/results/results_${sample}/latest.bad.reportTable.csv");


##################### Call script to summarise found phenotypes as XLS spreadsheet

print STDERR "\t +++ INFO: Combining Excel Tables\n";

#system("perl ${dir}scripts/GenomeChronicler_XLSX_fromTables.pl ${resultsdir}/results/results_${sample}/ ${resultsdir}/results/results_${sample}/${sample}_genotypes_${dtag}.xlsx");
system("python3 ${dir}scripts/GenomeChronicler_XLSX_fromTables.py ${resultsdir}/results/results_${sample}/ ${resultsdir}/results/results_${sample}/${sample}_genotypes_${dtag}.xlsx");


##################### Call LaTeX on the right template to produce the final report

print STDERR "\t +++ INFO: Compiling Genome Report\n";

system("cp $template ${resultsdir}/results/results_${sample}/${sample}_report_${dtag}.tex");
system("cp ${dir}templates/versionTable.txt ${resultsdir}/results/results_${sample}/");
system("cp ${dir}templates/GeneStructure.pdf ${resultsdir}/results/results_${sample}/");

my $TEMPLATETEX = "${sample}_report_".$dtag;

&runLatex();

sub runLatex {
    local $CWD = "${resultsdir}/results/results_${sample}";

	for(my $i = 0; $i < 3 ; $i++) {
		system("pdflatex -interaction=nonstopmode ${TEMPLATETEX} .tex 2> /dev/null >/dev/null");
	}
}

##################### Clean up temp files, both in the results folder and in the potentially extraneous tables produced for the report

print STDERR "\t +++ INFO: Cleaning up Temporary and Intermediate Files\n";

#note to self: maybe there will not be a BAM file if we are using a VCF file

system("rm -rf $BAM_file ${BAM_file}.bai");
system("rm -rf ${resultsdir}/results/results_${sample}/temp/");
system("rm -rf ${resultsdir}/results/results_${sample}/latest*.csv");
system("rm -rf ${resultsdir}/results/results_${sample}/versionTable.txt ${resultsdir}/results/results_${sample}/GeneStructure.pdf");
system("rm -rf ${resultsdir}/results/results_${sample}/${TEMPLATETEX}.out ${resultsdir}/results/results_${sample}/texput.log ${resultsdir}/results/results_${sample}/${TEMPLATETEX}.aux ${resultsdir}/results/results_${sample}/${TEMPLATETEX}.log ${resultsdir}/results/results_${sample}/${TEMPLATETEX}.tex");
system("rm -rf ${dir}GenomeChronicler_plot_generator_fromAncestry.Rout");

sleep (1);

$BAM_file =~ s/\.clean\.BAM//gi;
print STDERR "\n\t +++ DONE: Finished GenomeChronicler for file [ $BAM_file ] in ".( time() - $start_time )." seconds\n";






sub cleanBAMfile_noCHR() {

	my $filename = shift;

	#open(IN, "${dir}software/samtools view -H $filename |") or die "Could not open input file: $!\n";
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


	#system("${dir}software/samtools reheader $filename.tempHeader $filename > $filename.clean.BAM");
	#system("${dir}software/samtools index $filename.clean.BAM");	
	system("samtools reheader $filename.tempHeader $filename > $filename.clean.BAM");
	system("samtools index $filename.clean.BAM");

	#system("mv $filename.tempBAM $filename");
	system("rm $filename.tempHeader")
}


