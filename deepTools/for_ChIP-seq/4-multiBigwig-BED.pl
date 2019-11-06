#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;
## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "3-subtractInput/H3K4me3"
my $BEDdir_g = '';  ## such as "GenomicRegions/3_Enhancers/allPups"
my $output_g = '';  ## such as "5-multiBigwig-BED/H3K4me3"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Usage:
               perl  4-multiBigwig-BED.pl    [-version]    [-help]      [-in inputDir]     [-bed bedDir]     [-out outDir]
        For instance:
               perl  4-multiBigwig-BED.pl    -in 3-subtractInput/H3K4me3   -bed GenomicRegions/3_Enhancers/allPups    -out 5-multiBigwig-BED/H3K4me3    > 4-multiBigwig-BED.H3K4me3.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of input path that contains your BigWig files.  (no default)
        -bed bedDir         "bedDir"   is the name of input path that contains your BED files with genomic regions.  (no default)
        -out outDir         "outDir"   is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "version 0.1,  2019-11-05.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g    = '3-subtractInput/H3K4me3';     ## This is only an initialization value or suggesting value, not default value.
$BEDdir_g   = '3_Enhancers/allPups';         ## This is only an initialization value or suggesting value, not default value.
$output_g   = '5-multiBigwig-BED/H3K4me3';   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -bed   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  4-multiBigwig-BED.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-bed'     }   )     { $BEDdir_g = $args{'-bed'     }; }else{say   "\n -bed    is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$BEDdir_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input       Path:  $input_g
                Genomic  Regions:  $BEDdir_g
                Output      Path:  $output_g
        ###############################################################
\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say    "\n\n\n\n\n\n##################################################################################################";
say    "Running......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die; }
}

&myMakeDir($output_g);

opendir(my $DH_input_g, $input_g)   ||  die;
opendir(my $DH_bed_g,   $BEDdir_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
my @BED_Files_g  = readdir($DH_bed_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Checking all the necessary softwares in this step......" ;
sub printVersion  {
    my $software = $_[0];
    system("echo    '##############################################################################'  >> $output_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '#########$software'                                                              >> $output_g/VersionsOfSoftwares.txt   2>&1");
    system("$software                                                                                 >> $output_g/VersionsOfSoftwares.txt   2>&1");
    system("echo    '\n\n\n\n\n\n'                                                                    >> $output_g/VersionsOfSoftwares.txt   2>&1");
}
&printVersion("deeptools --version");
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting BigWig files in input folder ......";
my @BigWigfiles_g = ();
{
open(seqFiles_FH, ">", "$output_g/BigWig-Files.txt")  or  die; 
for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bw$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    next unless $inputFiles_g[$i] !~ m/^unpaired/;
    say    "\t......$inputFiles_g[$i]"; 
    $BigWigfiles_g[$#BigWigfiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBigWig file:  $inputFiles_g[$i]\n";
    say   seqFiles_FH  "BigWig file:  $inputFiles_g[$i]\n";
}

say   seqFiles_FH  "\n\n\n\n\n";  
say   seqFiles_FH  "All BigWig files:@BigWigfiles_g\n\n\n";
say        "\t\t\t\tAll BigWig files:@BigWigfiles_g\n\n";
my $num1 = $#BigWigfiles_g + 1;
say seqFiles_FH   "\nThere are $num1 BigWig files.\n";
say         "\t\t\t\tThere are $num1 BigWig files.\n";
}

my @BigWigfiles_g2 =  @BigWigfiles_g;    
for ( my $i=0; $i<=$#BigWigfiles_g2; $i++ ) { 
   $BigWigfiles_g2[$i] = "$input_g/$BigWigfiles_g2[$i]";   ## add path  
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Using multiBigwigSummary bin, plotCorrelation and plotPCA ......";

for ( my $i=0; $i<=$#BED_Files_g; $i++ ) { 
    next unless $BED_Files_g[$i] =~ m/\.bed$/;
    next unless $BED_Files_g[$i] !~ m/^[.]/;
    next unless $BED_Files_g[$i] !~ m/[~]$/;
    my $output_sub1 = "$output_g/$BED_Files_g[$i]";   
    &myMakeDir($output_sub1);
    system("multiBigwigSummary  BED-file    --bwfiles @BigWigfiles_g2    --BED $BEDdir_g/$BED_Files_g[$i]  --smartLabels   --numberOfProcessors max/2    --verbose    --outRawCounts $output_sub1/1-RawCounts.txt   --outFileName $output_sub1/1-results.npz    >> $output_sub1/1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/2-Correlation.heatmap.pearson.pdf       --outFileCorMatrix $output_sub1/2-Correlation.heatmap.pearson.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/3-Correlation.heatmap.spearman.pdf      --outFileCorMatrix $output_sub1/3-Correlation.heatmap.spearman.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/4-Correlation.scatterplot.pearson.pdf   --outFileCorMatrix $output_sub1/4-Correlation.scatterplot.pearson.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/1-results.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/5-Correlation.scatterplot.spearman.pdf  --outFileCorMatrix $output_sub1/5-Correlation.scatterplot.spearman.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz                                                        -o $output_sub1/6A-top1000-plotPCA.pdf                  --outFileNameData  $output_sub1/6A-top1000-plotPCA.txt                    --plotHeight 20   --plotWidth 20      >> $output_sub1/6A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz     --ntop 5000                                        -o $output_sub1/6B-top5000-plotPCA.pdf                  --outFileNameData  $output_sub1/6B-top10000-plotPCA.txt                   --plotHeight 20   --plotWidth 20      >> $output_sub1/6B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz     --ntop 0                                           -o $output_sub1/6C-topAll-plotPCA.pdf                   --outFileNameData  $output_sub1/6C-topAll-plotPCA.txt                     --plotHeight 20   --plotWidth 20      >> $output_sub1/6C-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz                     --rowCenter                        -o $output_sub1/7A-top1000-plotPCA.rowCenter.pdf        --outFileNameData  $output_sub1/7A-top1000-plotPCA.rowCenter.txt          --plotHeight 20   --plotWidth 20      >> $output_sub1/7A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz     --ntop 5000     --rowCenter                        -o $output_sub1/7B-top5000-plotPCA.rowCenter.pdf        --outFileNameData  $output_sub1/7B-top10000-plotPCA.rowCenter.txt         --plotHeight 20   --plotWidth 20      >> $output_sub1/7B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/1-results.npz     --ntop 0        --rowCenter                        -o $output_sub1/7C-topAll-plotPCA.rowCenter.pdf         --outFileNameData  $output_sub1/7C-topAll-plotPCA.rowCenter.txt           --plotHeight 20   --plotWidth 20      >> $output_sub1/7C-plotPCA-runLog.txt   2>&1");                               

}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END







