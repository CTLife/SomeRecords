#!/usr/bin/env perl
use strict;
use warnings;

## sudo apt-get install mysql-client-core-5.6
## ./fetchChromSizes.sh mm9 > mm9.chrom.sizes
## http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
## Convert bam to bigwig  for RNA-seq Data.


########## Keys and Values ##########
my %args = @ARGV;

########## Set Defaults ##########
my $input_g      = "";    ## Must be bam format
my $inputDir_g   = ""; 
my $output_g     = "";
my $outputDir_g  = "";
my $scale_g      = "";  ## number of reads


########## Get Arguments ##########
if ( exists $args{'-inputDir'}  ) { $inputDir_g  = $args{'-inputDir'};      }  
if ( exists $args{'-input'}     ) { $input_g     = $args{'-input'};         }   ## Must be bam format
if ( exists $args{'-outputDir'} ) { $outputDir_g = $args{'-outputDir'};     }
if ( exists $args{'-output'}    ) { $output_g    = $args{'-output'};        }
if ( exists $args{'-scale'}     ) { $scale_g     = $args{'-scale'};         }    ## number of reads


########### Conditions #############
$inputDir_g   =~ m/^\S+$/ or die;
$input_g      =~ m/^\S+$/ or die;   ## Must be bam format
$outputDir_g  =~ m/^\S+$/ or die;
$output_g     =~ m/^\S+$/ or die;
$scale_g      =~ m/^\S+$/ or die;   ## number of reads


######### Example ###########
# perl convertFormat.pl   -inputDir XXXXX   -input H2bGFP_week12_rep1.bam      -outputDir XXXXX     -output  H2bGFP_week12_rep1    -scale 64895321






my $bam      = "$inputDir_g/$input_g";
my $bedGraph = "$outputDir_g/$output_g.bedGraph";
my $bigwig   = "$outputDir_g/$output_g.bigwig";

$scale_g = 10**7/$scale_g;
print    "scale=$scale_g\n";
system("bedtools  genomecov   -bg   -split   -scale $scale_g      -ibam $bam    -g mm9.chrom.sizes   > $bedGraph");
sleep(2);
system("wigToBigWig  $bedGraph   mm9.chrom.sizes   $bigwig");

























