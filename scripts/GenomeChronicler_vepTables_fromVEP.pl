#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

#Sort out input
die "Usage: $0 VEPhtml OutputDir\n" if(@ARGV != 2);
my $filename = $ARGV[0];
my $outDir = $ARGV[1];



die "Major error here (File not found)\n" if(!-e $filename);

open(OUT, ">${outDir}/latest.summary.csv") or die "Could not open output file: $!\n";
print OUT "Feature\tCount\n";

open(IN, "grep \"gen_stats\" $filename |") or die "Could not open input file: $!\n";
my $recording = 0;
while (my $line = <IN>) {
    chomp($line);
    next if($line eq "");
    
    #Note for future generations: This is as silly as Perl gets, but it does the job (for now). I will re-write this with proper parsing at some point, or even better, use an off-the-shelve HTML parsing module.
    
    $line =~ s/<tr>/\n/g;
    $line =~ s/<\/table>/\n<\/table>\n/g;
    $line =~ s/<td>/\t/g;
    $line =~ s/<\/td>|<\/tr>//g;
    
    my @lines = split("\n",$line);
    for my $line (@lines) {
        if($line =~ m/stats_table/) {
            $recording += 1;
            next;
        }
        elsif($line =~ m/\<\/table\>/ and $recording == 2) {
            $recording = 0;
        }
        
        $line =~ s/^\s+//;
        $line =~ s/\s+$//;
        $line =~ s/\s{2,}/\t/g;
        
        next if($line =~ m/Variants processed/);
        next if($line =~ m/Lines of output written/);
        
        print OUT "$line\n" if($recording == 2);
    }
}

close IN;
close OUT;

#system("ln -sf ${filename}.summary.csv latest.summary.csv");


open(IN, "cat $filename | grep drawTable | grep \"\\[\" | ") or die "Could not open input file: $!\n";

while (my $line = <IN>) {
	chomp($line);
	next if($line eq "");

    my ($tabName,$tabString);
    
    if($line =~ m/drawTable\(\'(.+?)\'.+?\[(.+)\]/) {
        $tabName = $1;
        $tabString = $2;
    }
    
    my $counter = 0;
    my $summer = 0;
    my @table;
    while ($tabString =~ m/\[(.+?)\]/g) {
        my $pair = $1;
        $pair =~ s/[']//g;
        $pair =~ s/\_/ /g;
        
        $pair =~ s/\ //g if($counter == 0); #Remove spaces in title only
        my @split = split(",",$pair);
        push(@table,\@split);
        
        $summer += $split[1] if($counter >0);
        $counter++;
    }
    
    open(OUT, ">${outDir}/latest.${tabName}.csv") or die "Could not open output file: $!\n";
    
    my $head = shift(@table);
    $head->[0] = "Label";
    print OUT join(",",@{$head})."\n";
    
    my $other = 0;
    for my $d (@table) {
        my $perc = sprintf("%0.2f",(($d->[1]/$summer)*100));
        if($perc > 1) {
            $d->[0] = ucfirst($d->[0]);
            $d->[1] = $perc;
            print OUT join(",",@{$d})."\n";
        }
        else {
            $other += $d->[1];
        }
    }
    if($other > 0) {
        my $perc = sprintf("%0.2f",(($other/$summer)*100));
        print OUT "Others,$perc\n";
    }
    close OUT;
    
    #system("ln -sf ${filename}.${tabName}.csv latest.${tabName}.csv");
    
    #die "$tabName\n$tabString\n";
}

close IN;
