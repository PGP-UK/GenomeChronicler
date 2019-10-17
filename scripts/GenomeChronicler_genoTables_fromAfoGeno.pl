#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI;

my $version = "19-257";

# Note: It is quite clear this code can be optimised, but it appears to work as it is:

my %tmpBlacklist = ("rs10156191"=>1,"rs12094543"=>1,"rs1667255"=>1,"rs17672135"=>1,"rs6445975"=>1,"rs7659604"=>1);
    

print STDERR "This is script [$0] part of the report generator pipeline version [$version];\n";

#Sort out input
die "Usage: $0 MassagedGenotypeVCF OutputDirectory\n" if(@ARGV != 2);
my $filename = $ARGV[0];
#my $sampleFR = $ARGV[1];
my $outdir = "./";
$outdir = $ARGV[1] if(defined($ARGV[1]));

###print STDERR "ToDo:Backport the extra script back here to handle Lucia\'s new VCF files including RSID.\n\n";


my $dir="/GenomeChronicler/";

my $driver   = "SQLite";

my $database = "${dir}reference/snpedia.db";
my $dsn = "DBI:$driver:dbname=$database";
my $dbh = DBI->connect($dsn, "", "", { RaiseError => 1 }) or die $DBI::errstr;
print "Opened database successfully [SNPedia]\n";
#my $stmt = qq(select * from data where id=\"?\");
#my $sth = $dbh->prepare( $stmt );
#my $sth1a = $dbh->prepare( 'select * from data join strand where data.id=? and strand.id=data.id' );
my $sth1c = $dbh->prepare( 'select * from data where data.id=?' );
my $sth1f = $dbh->prepare( 'select * from flagged where id=?' );
my $sth1s = $dbh->prepare( 'select * from strand where id=?' );
my $sth1g = $dbh->prepare( 'select * from genoset where id=?' );

#my $database2 = "results_${sampleFR}/${sampleFR}.db";
#my $dsn2 = "DBI:$driver:dbname=$database2";
#my $dbh2 = DBI->connect($dsn2, "", "", { RaiseError => 1 })
#or die $DBI::errstr;
#print "Opened database successfully [VolunteerVCF]\n";
#my $sth2 = $dbh2->prepare('select * from data where chr=? and coord=?');


my $database3 = "${dir}reference/gnomad.db";
my $dsn3 = "DBI:$driver:dbname=$database3";
my $dbh3 = DBI->connect($dsn3, "", "", { RaiseError => 1 })
or die $DBI::errstr;
print "Opened database successfully [GnomAD]\n";
my $sth3 = $dbh3->prepare('select * from data where rsid=?');


my $database4 = "${dir}reference/getevidence.db";
my $dsn4 = "DBI:$driver:dbname=$database4";
my $dbh4 = DBI->connect($dsn4, "", "", { RaiseError => 1 })
or die $DBI::errstr;
print "Opened database successfully [GetEvidence]\n";
my $sth4 = $dbh4->prepare('select * from data where dbsnp_id=?');


my $database5 = "${dir}reference/clinvar.db";
my $dsn5 = "DBI:$driver:dbname=$database5";
my $dbh5 = DBI->connect($dsn5, "", "", { RaiseError => 1 })
or die $DBI::errstr;
print "Opened database successfully [ClinVar]\n";
my $sth5 = $dbh5->prepare('select * from data where rsid=?');



#Input Genoset Data


#Note: This is probably useless and never worked 100%, so it is ready to be deleted (just create the array from sorted keys from the hash below).
open(IN, "${dir}reference/genosetDependencies.txt") or die "Could not open input file: $!\n";
my @allGenosets;
while (my $line = <IN>) {
    chomp($line);
    next if($line eq "");
    my @data = split("\t", $line);
    
    push(@allGenosets,$data[0]);
}
close IN;
#die Dumper(\@allGenosets);


open(IN, "${dir}reference/parsedGenosets.txt") or die "Could not open input file: $!\n";
my %genosets;
my %genotypes;
while (my $line = <IN>) {
    chomp($line);
    next if($line eq "");
    my @data = split("\t", $line);
    
    my $name = lc(shift(@data));
    $genosets{$name} = \@data;
    #die Dumper(\@data);
}
close IN;

#<STDIN>;
#die Dumper(\%genosets);


open(IN, "$filename") or die "Could not open input file [114] [[$filename]]: $!\n";

open(GOOD, ">${outdir}/latest.good.reportTable.csv") or die "Could not open output file: $!\n";
#print GOOD "Effect,Identifier,Magnitude,Genotype,Summary,GnomAD,GetEvidence,ClinVar\n";
print GOOD "Mag.,Identifier,Genotype,Summary,GnomAD,GetEvidence,ClinVar\n";
close GOOD;
open(GOOD, " | sort | uniq | sort -t \',\' -k1,1nr >> ${outdir}/latest.good.reportTable.csv") or die "Could not open output file: $!\n";

open(BAD, ">${outdir}/latest.bad.reportTable.csv") or die "Could not open output file: $!\n";
print BAD "Mag.,Identifier,Genotype,Summary,GnomAD,GetEvidence,ClinVar\n";
close BAD;
open(BAD, " | sort | uniq | sort -t \',\' -k1,1nr >> ${outdir}/latest.bad.reportTable.csv") or die "Could not open output file: $!\n";



while (my $line = <IN>) {
	chomp($line);
	next if($line eq "");
    next if($line =~ m/^#/);
        
        
    my $debugBuffer = "";
	
    #OLD WAS: 1_54713_-/TTTC 1:54712-54713 TTTC ENSG00000268020 ENST00000606857 Transcript downstream_gene_variant - - - - - rs568927205 IMPACT=MODIFIER;DISTANCE=1400;STRAND=1;SYMBOL=OR4G4P;SYMBOL_SOURCE=HGNC;HGNC_ID=14822;BIOTYPE=unprocessed_pseudogene;CANONICAL=YES;AFR_MAF=TTTC:0.5083;AMR_MAF=TTTC:0.5677;EAS_MAF=TTTC:0.6776;EUR_MAF=TTTC:0.6193;SAS_MAF=TTTC:0.544

    #NEW IS: 7 92383887 92383888 rs10 7 92383888 A C 1/1 C C

    
    #Even newer... (NOTE: Need to adapt to this).
    #MT	16519	rs3937033	T	C	103537	.	AC=2;AF=1;AN=2;BaseQRankSum=-0.543;ClippingRankSum=0.19;DP=3147;FS=0;MLEAC=2;MLEAF=1;MQ=59.94;MQRankSum=-0.622;QD=33.24;ReadPosRankSum=-0.088;SOR=0.728	GT:AD:DP:GQ:PL	1/1:5,3110:3118:99:103565,9151,0
    
    my @data = split("\t", $line);
    $data[11] = "" if(!defined($data[11]));
    
    $debugBuffer .= "\n+++\n@data\n";
    
    #die "@data\n";
    
    
    #This is not needed
    #my $chr = $data[0];
    #my $coord = $data[5];
    
    my $strand = "plus";

    my $snp = $data[3];

    my $rvS = $sth1s->execute($snp) or die $DBI::errstr;
    die $DBI::errstr if($rvS < 0);
    
    my $counter = 0;
    while(my @row = $sth1s->fetchrow_array()) {
        $strand = $row[1];
        $counter++;
    }
    die "Oh my Gawd!!! DNA has two strands (and it isn't a good thing in this case)" if($counter > 1);
    
    
    
    
    my $rvFlag = $sth1f->execute($snp) or die $DBI::errstr;
    die $DBI::errstr if($rvFlag < 0);
    my $flagID = 0;
    while(my @row = $sth1f->fetchrow_array()) {
        $flagID = 1;
        #die "Found a flag! Hurray! Flag! I found a Flag! Will someone listen to me? I have found a flag!!! [$snp]\n";
    }
    
    #die "Found my candidate SNP" if($snp eq "rs141935559");

    
    #die "Do something here about the flag if it is raised and genotype = 0/0, probably just next it... = [geno] $data[8] = [extra] $data[11]";
    
    if(($flagID == 1) and (($data[8] eq "0/0") or ($data[8] eq "./."))) {
        #Uncomment below for DEBUG: Indels that are flagged and don't have info or are reference.
        #print STDERR "POTENTIAL ERROR: Ups, found a bugger. Kill the bugger... next! [$snp]\n";
        next;
    }
    
    $data[9]  =~ tr/ACTGactg/TGACtgac/ if(lc($strand) eq "minus");
    $data[10] =~ tr/ACTGactg/TGACtgac/ if(lc($strand) eq "minus");
    my ($genotype,$extra) =  ([$data[9],$data[10]],$data[11]);
    
    $genotypes{$snp}{$extra} = 1 if($extra ne "");
    $genotypes{$snp}{$genotype->[0]} = 1;
    $genotypes{$snp}{$genotype->[1]} = 1;
    $genotypes{$snp}{"$genotype->[0];$genotype->[1]"} = 1;
    
    
    
    my $rv = $sth1c->execute($snp) or die $DBI::errstr;
    if($rv < 0){
        print $DBI::errstr;
    }
    my %alleles;
    my $countGenotype = 0;
    
    
    while(my @row = $sth1c->fetchrow_array()) {

        $debugBuffer .= "ID = ". $row[0] . "\n";
        $debugBuffer .= "REPUTE = ". $row[1] ."\n";
        $debugBuffer .= "MAGNITUDE = ". $row[2] ."\n";
        $debugBuffer .= "ALLELE1 =  ". $row[3] ."\n";
        $debugBuffer .= "ALLELE2 =  ". $row[4] ."\n";
        $debugBuffer .= "SUMMARY =  ". $row[5] ."\n";
        $debugBuffer .= "RSID =  ". $row[6] ."\n";
        $debugBuffer .= "IID =  ". $row[7] ."\n";
        #$debugBuffer .= "ID2 =  ". $row[8] ."\n";
        #$debugBuffer .= "STRAND =  ". $row[9] ."\n";

        
        #die "DEBUUUG_cycle: $line\n $debugBuffer\n" if($line =~ m/rs3928306/);

        #$strand = $row[9];
        
        #Doing a full summary string cleanup (now in a neat subroutine to be re-used in genosets).
        $row[5] = cleanSummaryString($row[5]);
        
        my $a1 = $row[3];
        my $a2 = $row[4];
        
        #I'm swapping the genotypes now, and not the alleles here.
        #$a1 =~ tr/ACTGactg/TGACtgac/ if(lc($strand) eq "minus");
        #$a2 =~ tr/ACTGactg/TGACtgac/ if(lc($strand) eq "minus");
        
        next if(!defined($row[2]) or ($row[2] eq 0)); #Excluding magnitude 0
     
        $alleles{$a1}{$a2} = [$row[1],$row[0],$row[2],"($row[3];$row[4])",$row[5]];
        $alleles{$a2}{$a1} = $alleles{$a1}{$a2};
        
        $countGenotype++;
    }
    


    
        if($countGenotype) {
            
            if(defined($alleles{$genotype->[0]}{$genotype->[1]})) {
                
            
                #Note: This is a bit silly but it is happening here instead of pre-filtering the database so that the SNPs can be used for the genoset calculation, even if we are not going to report them to the user. These can also be filtered even later in the process.                 
                
                #Implement some filters here
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^\s*$/);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Common/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^None/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Normal/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Average/);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Benign most likely/);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Ancestral value/);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Unaffected Genotype/i);
                
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Benign variant/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^L22- S142-/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Depends on/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Extensive metabolizer/i);

                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Typical BuChE/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Complex; generally normal risk/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/common in clinvar/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Major allele, normal risk/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Slight if any/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Major allele, normal risk/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Likely to be benign/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Most common genotype/i);


                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Most likely a benign polymorphism/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Benign \(harmless\) variant/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^1.3x risk$/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/normal risk of migraine/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/More likely to go bald; common/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Most likely benign polymorphism/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Slightly increased lifespan?/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/^Benign polymorphism/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/No increased risk of /i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Likely to be a benign variant/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Carrier of a benign change/i);
                next if($alleles{$genotype->[0]}{$genotype->[1]}->[4] =~ m/Classified as benign in ClinVar/i);


                               
                
                my $rsid = $alleles{$genotype->[0]}{$genotype->[1]}->[1];
                my $GnomAD = generate_GnomAD_url($rsid);
                my $getevidence = generate_GetEvidence_url($rsid);
                my $snpedia = generate_SNPedia_url($rsid);
                my $clinvar = generate_ClinVar_url($rsid);
                                
                next if(($GnomAD eq "") and ($getevidence eq "") and ($clinvar eq "")); #No links, no joy
                
                #Quick temporary hack (that is still here)
                next if(defined($tmpBlacklist{$rsid}));
                
                #Quick hack here to replace the RSid by a SNPedia link
                $alleles{$genotype->[0]}{$genotype->[1]}->[1] = $snpedia;
                
                #Generate link for summary
                $alleles{$genotype->[0]}{$genotype->[1]}->[4] = generate_Summary_url($rsid,$alleles{$genotype->[0]}{$genotype->[1]}->[4]);

                
                #my $tmpString = join(",",@{$alleles{$genotype->[0]}{$genotype->[1]}});
                my $tmpString = $alleles{$genotype->[0]}{$genotype->[1]}->[2].",".$alleles{$genotype->[0]}{$genotype->[1]}->[1].",".$alleles{$genotype->[0]}{$genotype->[1]}->[3].",".$alleles{$genotype->[0]}{$genotype->[1]}->[4];
                
                #my $tmpString = join("\t",@{$alleles{$genotype->[0]}{$genotype->[1]}});
                $tmpString =~ s/[^[:ascii:]]//g; #Clearing long text from unprintable characters caused by codepage mismatch along the way... try to find a better way to do this.
                
                next if(!defined($alleles{$genotype->[0]}{$genotype->[1]}->[2]) or $alleles{$genotype->[0]}{$genotype->[1]}->[2] eq "");
                
                if(!defined($alleles{$genotype->[0]}{$genotype->[1]}->[0]) or $alleles{$genotype->[0]}{$genotype->[1]}->[0] eq "") {
                    
                    #Uncomment for DEBUG
                    #print STDERR "\nWARNING: Go here and with a fair sense of justice determine is this is good or bad\n$tmpString\n\n";
                    next;
                }
                elsif($alleles{$genotype->[0]}{$genotype->[1]}->[0] eq "Good") {
                    print GOOD "$tmpString,$GnomAD,$getevidence,$clinvar\n";
                }
                elsif($alleles{$genotype->[0]}{$genotype->[1]}->[0] eq "Bad") {
                    print BAD "$tmpString,$GnomAD,$getevidence,$clinvar\n";
                }
                
                

            }
            elsif ($countGenotype > 2) {
                print STDERR "IMPORTANT: $debugBuffer\n";
                print STDERR "IMPORTANT: $debugBuffer\n";
                print STDERR "IMPORTANT: Please debug this here as I couldn't find the right alleles [$strand] [$countGenotype]\n\n\n";
                
            }
            
        }
#        else { #Uncomment this and run to debug the variants that are being discarded because their rsIDs are not in SNPedia
#            for my $lines (@data) {
#                print STDERR "WARNING --- NOT FOUND --- \t$lines\n";
#            }
#            #die;
#        }
    
}

close IN;

close GOOD;
close BAD;



##Hacking gs269...
#$genotypes{"rs7412"}{"C;T"}=1;
#$genotypes{"rs429358"}{"T;T"}=1;
#die Dumper(\%genotypes);
#die Dumper($genotypes{"rs3928306"});



############################### Having a go at processing the genosets...

my %positiveGenosets;
for my $genoset (@allGenosets) {
    
    #next if($genoset ne "gs269");
    #next if($genoset ne "gs237");
    my $result = processGenoset($genoset);
    
}

#Slight debug here... Revisit this after new genotypes used as input.
#die Dumper(\%genotypes);
print STDERR "\t +++ POSITIVES:\n".Dumper(\%positiveGenosets);

open(GENO, ">${outdir}/latest.genoset.reportTable.csv") or die "Could not open output file: $!\n";
print GENO "Magnitude,Identifier,Summary\n";
close GENO;

open(GENO, " | sort -t \",\" -k 1,1nr >>${outdir}/latest.genoset.reportTable.csv") or die "Could not open output file: $!\n";

for my $geno (sort(keys(%positiveGenosets))) {
    
    #print "Trying to get $geno\n\n";
    my $rvG = $sth1g->execute($geno) or die $DBI::errstr;
    if($rvG < 0){
        die "GENOSET DB ERROR: ".$DBI::errstr;
    }
    
    while(my @row = $sth1g->fetchrow_array()) {
        #print STDERR Dumper(\@row);
        if($row[2] ne "0") {
            print GENO "$row[2],".generate_SNPedia_url($row[0]).",".generate_Summary_url($row[0],cleanSummaryString($row[3]))."\n";
        }
    }
    
}

close GENO;



#system("ln -sf ${filename}.genoset.reportTable.csv latest.genoset.report.csv");
#system("ln -sf ${filename}.good.reportTable.csv latest.good.report.csv");
#system("ln -sf ${filename}.bad.reportTable.csv latest.bad.report.csv");


print "Operation done successfully\n";
$dbh->disconnect();










sub processGenoset {
    my $genoset = shift;
    
    
    #print STDERR Dumper($genosets{$genoset})."\n";
    
    if(!defined($genosets{$genoset}->[0])) {
        print STDERR "Logic conditional not found for genoset $genoset, please double-check this genoset has been deprecated\n";
        return "0";
    }
    
    my $logic = $genosets{$genoset}->[0];
    my @vars = split(",",$genosets{$genoset}->[1]);
    
    my @v = ();
    my @g = ();
    
    
    #BIG FAT WARNING: DO SOMETHING CLEVER ABOUT THE iIDs OTHERWISE WE WILL BE IN TROUBLE
    
    #Fill in vars
    for my $var (@vars) {
        #print STDERR "\tFilling $var\n";
        my @varComponents = split("=",$var);
        
        eval("$varComponents[0] = \"".evaluateGenotype($varComponents[1])."\";")
        
    }
    
    #Evaluate Genoset
    my $result = 0;
    #eval("print STDERR \"$logic\n\"");
    eval("\$result = 1 if($logic);");
    
    #print STDERR "RESULT HERE [$genoset] IS: $result for [$logic] and [[@v <-v  g-> @g]]\n";

    #Debug break
    #die Dumper(\@v,\@g,$result) if($result eq "1");
    
    #Put data back into genotypes for future use
    $genotypes{$genoset} = $result;
    
    #If positive, add to positiveGenosets (note: I know I shouldn't change global variables from within subroutines... I'll pay for my sins someday.)
    $positiveGenosets{$genoset} = 1 if ($result eq "1");

    return $result;
}




sub evaluateGenotype {
    my $query = shift;
    
    
    my $var;
    my $value;
    my $result = "NA";
    
    
    if($query =~ m/(.+)\((.+)\)/) {
        my $snp = $1;
        my $geno = $2;
        
        #die "DEBUG QUERY [362]: $query [$snp] [$geno]\n" if($query =~ m/rs3928306/);

        #print STDERR "Checking for $snp $geno".Dumper($genotypes{$snp})."\n\n";
        
        if(!defined($genotypes{$snp})) {
            print STDERR "Some big problem here hum? [$snp]\n";
        }
        
        
        
        if(defined($genotypes{$snp}{$geno})) {
            return 1;
        }
        else {
            return 0;
        }
    }
    elsif($query =~ m/(gs.+)/i) {
        my $genoset = $1;
        
        #print STDERR "Checking for genoset: $genoset".Dumper($genotypes{$genoset})."\n\n";

        if(defined($genotypes{$genoset})) {
            return $genotypes{$genoset};
        }
        else {
            #return 1; #Note: this is stupid and false...
            #die "Dependency fail, genoset $genoset was needed and not found\n";
            
            #print STDERR "WARNING: Genoset $genoset should have been previously computed... (no harm done unless there is a dependency cycle, at which point you have a much bigger problem that this to solve)\n";
            return processGenoset($genoset); #Sorry for the silly recursiveness... I hope there are no cycles and otherwise I'm in trouble...
        }
    }
    else {
        die "I'm sorry Dave, I'm afraid I can't do that [$query]\n";
    }
    
    
    return $result;
}






#This was updated on 17 of June 2016 in preparation for the input format change. At the moment this function is not called at all.
sub buildGenotype {
    my @alleles = ($_[0],split(",",$_[1]));
    my @extras = ("");
    my $gen = $_[2];
    
    
    my $l1 = length($alleles[0]);
    my $extra = "";
    
    
    for (my $i = 1; $i < @alleles; $i++ ) {
        
        $extras[$i] = "";
        
        my $l2 = length($alleles[$i]);
        if($l1 > 1 or $l2 > 1) {
            
            ##Debug statement
            #print STDERR "There is an indel here [@alleles] [$gen]\n";
            
            
            if($l2 > $l1) {
                $alleles[0] = "-";
                $alleles[$i] =~ s/^.//s ;
                #die "There is an insertion here [@alleles] [$gen]\n";
                $extras[$i] = "I";
            }
            else {
                $alleles[0] =~ s/^.//s ;
                $alleles[$i] = "-";
                #die "Implement deletions here please  [@alleles] [$gen]\n";
                $extras[$i] = "D";
            }
        }
    }
    
    
    #    die "Genotype debug [$gen]\n" if($gen ne "0/1" and $gen ne "1/1"); #This is to account for odd genotypes. None have been encountered so far but they exist... #Now solved in this version I think
    
    my @res = map{$alleles[$_]} split('/',$gen);
    $extra = join(";",map{$extras[$_]} split('/',$gen));
    $extra =~ s/^;//;
    $extra =~ s/;$//;
    
    #print STDERR "@alleles\t$gen\t@res\n";
    
    #die "Just come back here to buildGenotype and check this case is working as it should\n=====\nINPUT:@_ -- @alleles -- $gen == @res --- $extra\n=====\n" if(length($extra) > 1);
    
    print STDERR "\n=====\nINPUT:@_ -- @alleles -- $gen == @res\n=====\n" if($_[2] =~ m/2/);
    
    return \@res,$extra;
}



#sub buildGenotype.OLD {
#    my @alleles = ($_[0],split(",",$_[1]));
#    my @extras = ("");
#    my $gen = $_[2];
#    
#    
#    my $l1 = length($alleles[0]);
#    my $extra = "";
#    
#
#    for (my $i = 1; $i < @alleles; $i++ ) {
#        
#        $extras[$i] = "";
#        
#        my $l2 = length($alleles[$i]);
#        if($l1 > 1 or $l2 > 1) {
#            
#            ##Debug statement
#            #print STDERR "There is an indel here [@alleles] [$gen]\n";
#            
#            
#            if($l2 > $l1) {
#                $alleles[0] = "-";
#                $alleles[$i] =~ s/^.//s ;
#                #die "There is an insertion here [@alleles] [$gen]\n";
#                $extras[$i] = "I";
#            }
#            else {
#                $alleles[0] =~ s/^.//s ;
#                $alleles[$i] = "-";
#                #die "Implement deletions here please  [@alleles] [$gen]\n";
#                $extras[$i] = "D";
#            }
#        }
#    }
#
#    
#    #    die "Genotype debug [$gen]\n" if($gen ne "0/1" and $gen ne "1/1"); #This is to account for odd genotypes. None have been encountered so far but they exist... #Now solved in this version I think
#    
#    my @res = map{$alleles[$_]} split('/',$gen);
#    $extra = join(";",map{$extras[$_]} split('/',$gen));
#    $extra =~ s/^;//;
#    $extra =~ s/;$//;
#    
#    #print STDERR "@alleles\t$gen\t@res\n";
#    
#    
#    die "Just come back here to buildGenotype and check this case is working as it should\n=====\nINPUT:@_ -- @alleles -- $gen == @res --- $extra\n=====\n" if(length($extra) > 1);
#    
#    #print STDERR "\n=====\nINPUT:@_ -- @alleles -- $gen == @res\n=====\n" if($_[2] =~ m/2/);
#
#    
#    return \@res,$extra;
#}
#


sub cleanSummaryString {
    
    my $row = shift;
    
    #Removing commas from summary so we can use a CSV format for the LaTeX table later on.
    $row =~ s/\,/\:/g;
    $row =~ s/\&/and/g;
    
    #Tidying up the summaries to appease my OCD
    $row = ucfirst($row);
    
    #Sorry, but I'll have to cut on the info here...
    $row = substr($row,0,47);
    $row .= "..." if(length($row) == 47);
    
    
    return $row;
}




sub generate_GnomAD_url {
    my $rsid = shift;
    
    my $rv3 = $sth3->execute($rsid) or die $DBI::errstr;
    if($rv3 < 0){
        die $DBI::errstr;
    }

    my $counter = 0;
    while(my @row = $sth3->fetchrow_array()) {
        $counter++;
    }

    if($counter) {
        #return "http://exac.broadinstitute.org/dbsnp/$rsid";
        return '\href{http://gnomad.broadinstitute.org/awesome?query='.$rsid.'}{Link}';
    }
    else {
        return "";
    }
    
}


sub generate_GetEvidence_url {
    my $rsid = shift;
    
    my $rv4 = $sth4->execute($rsid) or die $DBI::errstr;
    if($rv4 < 0){
        die $DBI::errstr;
    }
    
    my $counter = 0;
    while(my @row = $sth4->fetchrow_array()) {
        $counter++;
    }
    
    if($counter) {
        #return "http://evidence.pgp-hms.org/$rsid";
        return "\\href{http://evidence.pgp-hms.org/$rsid}{Link}";
    }
    else {
        return "";
    }
    
}

sub generate_SNPedia_url {
    my $rsid = shift;
    
    #return "$rsid";
    #return "https://www.snpedia.com/index.php/".ucfirst($rsid);
    return '\href{https://www.snpedia.com/index.php/'.ucfirst($rsid).'}{'.$rsid.'}';
}




sub generate_Summary_url {
    my $rsid = shift;
    my $summary = shift;
    
    return '\href{https://www.snpedia.com/index.php/'.ucfirst($rsid).'}{'.$summary.'}';
}


sub generate_ClinVar_url {
    my $rsid = shift;
    
    my $rv5 = $sth5->execute($rsid) or die $DBI::errstr;
    if($rv5 < 0){
        die $DBI::errstr;
    }

    my $counter = 0;
    while(my @row = $sth5->fetchrow_array()) {
        $counter++;
    }

    if($counter) {
        #return "http://www.ncbi.nlm.nih.gov/clinvar/?term=$rsid";
        return "\\href{http://www.ncbi.nlm.nih.gov/clinvar/?term=$rsid}{Link}";
    }
    else {
        return "";
    }

}



