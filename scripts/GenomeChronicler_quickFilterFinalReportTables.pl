#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use File::Basename;

#my @files = qw/latest.good.report.csv latest.bad.report.csv/;

my $filename = $ARGV[0];

my $DEBUG = 0;

#Sort out input
my $id = time();

#for my $filename (@files) {

	open(IN, "$filename") or die "Could not open input file: $!\n";
	open(OUT, ">$filename.temp.$id.txt") or die "Could not open output file: $!\n";

	my $head = <IN>;
	print OUT "$head";

	while (my $line = <IN>) {
		chomp($line);
		next if($line eq "");
		my @data = split(",", $line);

		print Dumper(\@data)."" if($DEBUG == 1);

		if($data[0] == 0) {
			print "WARNING: Magnitude 0\n" if($DEBUG == 1);
			next;
		}

		if(($data[1] !~ m/\{rs/) or ($data[1] =~ m/\{i/)) {
			print "WARNING: Non RS identifier\n" if($DEBUG == 1);
			next;
		}

		if(scalar(@data)<=4){
			print STDERR "WARNING: No links found\n" if($DEBUG == 1);
			next;
		}

		print "\n" if($DEBUG == 1);

		print OUT $line."\n";
#		die Dumper(\@data)."\n";

	}

	close IN;
	close OUT;

	system("cat $filename.temp.$id.txt > $filename");
	system("rm -rf $filename.temp.$id.txt");
#}
