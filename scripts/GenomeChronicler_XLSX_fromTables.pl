#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Excel::Writer::XLSX;

#Sort out input
print STDERR "Usage: $0 BaseDirectory OutFilename\n" if(!@ARGV);
my $basefilename = $ARGV[0];
my $outfilename = $ARGV[1];

my $workbook = Excel::Writer::XLSX->new($outfilename);

$workbook->set_properties(
title    => 'Excel format Genome Report',
author   => 'Jose Afonso Guerra-Assuncao',
comments => 'Generated for the Personal Genomes Project - United Kingdom Study',
);


my %extentions = ("Possibly Beneficial"=>"latest.good.reportTable.csv","Possibly Harmful"=>"latest.bad.reportTable.csv","Genosets"=>"latest.genoset.reportTable.csv");

my $urlformat = $workbook->add_format( color => 'black', underline => 1 );
my $headerformat = $workbook->add_format( color => 'black', bold => 1 );



for my $ext (sort(keys(%extentions))) {

    my $worksheet = $workbook->add_worksheet($ext);
    my $cr = 0;
    
    my @maxSize;
    
    my $filename = $basefilename.$extentions{$ext};
    #print STDERR "Adding $filename to XLSX\n";
    open(IN, "$filename") or die "Could not open input file: $!\n";

    while (my $line = <IN>) {
        chomp($line);
        next if($line eq "");
        my @data = split(",", $line);

        for (my $cc = 0; $cc < @data; $cc++) {
            my $strLength = 0;
            
            if($data[$cc] =~ m/href\{(.+)\}\{(.+)\}/) {
                #print STDERR "Something $data[$cc]\n";
                $worksheet->write_url($cr, $cc, $1, $urlformat, $2);
                $strLength = length($2);
            }
            else {
                if($cr == 0) {
                    $worksheet->write($cr, $cc, $data[$cc], $headerformat);
                }
                else {
                    $worksheet->write($cr, $cc, $data[$cc]);
                }
                $strLength = length($data[$cc]);
            }

        
            if(defined($maxSize[$cc])) {
                if($strLength > $maxSize[$cc]) {
                    $maxSize[$cc] = $strLength;
                }
            }
            else {
                $maxSize[$cc] = $strLength;
            }

        }
        
        $cr++;
    }

    for (my $cc = 0; $cc < @maxSize; $cc++) {

        $worksheet->set_column($cc, $cc, $maxSize[$cc]);
        
    }
    
    close IN;
}

$workbook->close() or die "Error closing file: $!";

