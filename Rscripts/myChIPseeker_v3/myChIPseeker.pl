#!/usr/bin/env perl
use strict;
use warnings;



######################################################################################################################
## Keys and Values
my %args = @ARGV;

## Initialize  Variables
my $inDir1  = '5-peaksOverlap-pups-rename';          
my $inDir2  = 'ART-vs-NC-H3K4me1';   
my $outDir  = "6-overlapPupsPeaks-annotation";

## Get Arguments
$inDir1  = $args{'-inDir1'  };          
$inDir2  = $args{'-inDir2'  };   
$outDir  = $args{'-outDir'  };

## Example
#  perl  myChIPseeker.pl    -inDir1 5-peaksOverlap-pups-rename     -inDir2  ART-vs-NC-H3K4me1       -outDir 6-overlapPupsPeaks-annotation    
######################################################################################################################



opendir(DIRHANDLE, "$inDir1/$inDir2")  or die; 
system("mkdir -p  $outDir");



while ( my $file1 = readdir(DIRHANDLE) ) {       
        next unless $file1 !~ m/^[.]/; 
        next unless $file1 !~ m/[~]$/;
        next unless $file1 =~ m/\.bed$/;
        system("mkdir -p  $outDir/$inDir2");
        system("Rscript  myChIPseeker.hg38.R     $inDir1     $inDir2/$file1    $outDir    > $outDir/$inDir2/$file1.runLog  2>&1 ");                   
} 
























