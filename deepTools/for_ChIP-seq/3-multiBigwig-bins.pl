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
my $output_g = '';  ## such as "4-multiBigwig-bins/H3K4me3"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Usage:
               perl  3-multiBigwig-bins.pl    [-version]    [-help]      [-in inputDir]    [-out outDir]
        For instance:
               perl  3-multiBigwig-bins.pl    -in 3-subtractInput/H3K4me3   -out 4-multiBigwig-bins/H3K4me3    > 3-multiBigwig-bins.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of input path that contains your BigWig files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "version 0.1,  2019-10-31.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '3-subtractInput/H3K4me3';     ## This is only an initialization value or suggesting value, not default value.
$output_g = '4-multiBigwig-bins/H3K4me3';  ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  3-multiBigwig-bins.pl  -help' \n";
    exit 0;
}

## Get Arguments
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input       Path:  $input_g
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

opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
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

{
     my $output_sub1 = "$output_g/multiBigwigSummary_bin/10000bp";   
     &myMakeDir($output_sub1);
    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 10000   --numberOfProcessors max/2    --verbose    --outRawCounts $output_sub1/A.1-RawCounts.10000bp-Bin.txt   --outFileName $output_sub1/A.1-results.10000bp.npz    >> $output_sub1/A.1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/A.1-results.10000bp.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/A.2-Correlation.heatmap.pearson.10000bp.pdf       --outFileCorMatrix $output_sub1/A.2-Correlation.heatmap.pearson.10000bp.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/A.2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/A.1-results.10000bp.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/A.3-Correlation.heatmap.spearman.10000bp.pdf      --outFileCorMatrix $output_sub1/A.3-Correlation.heatmap.spearman.10000bp.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/A.3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/A.1-results.10000bp.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/A.4-Correlation.scatterplot.pearson.10000bp.pdf   --outFileCorMatrix $output_sub1/A.4-Correlation.scatterplot.pearson.10000bp.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/A.4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/A.1-results.10000bp.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/A.5-Correlation.scatterplot.spearman.10000bp.pdf  --outFileCorMatrix $output_sub1/A.5-Correlation.scatterplot.spearman.10000bp.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/A.5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/A.1-results.10000bp.npz                                                        -o $output_sub1/A.6A-top1000-plotPCA.10000bp.pdf                  --outFileNameData  $output_sub1/A.6A-top1000-plotPCA.10000bp.txt                    --plotHeight 20   --plotWidth 20      >> $output_sub1/A.6A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/A.1-results.10000bp.npz     --ntop 10000                                       -o $output_sub1/A.6B-top10000-plotPCA.10000bp.pdf                 --outFileNameData  $output_sub1/A.6B-top10000-plotPCA.10000bp.txt                   --plotHeight 20   --plotWidth 20      >> $output_sub1/A.6B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/A.1-results.10000bp.npz     --ntop 0                                           -o $output_sub1/A.6C-topAll-plotPCA.10000bp.pdf                   --outFileNameData  $output_sub1/A.6C-topAll-plotPCA.10000bp.txt                     --plotHeight 20   --plotWidth 20      >> $output_sub1/A.6C-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/A.1-results.10000bp.npz                     --rowCenter                        -o $output_sub1/A.7A-top1000-plotPCA.10000bp.rowCenter.pdf        --outFileNameData  $output_sub1/A.7A-top1000-plotPCA.10000bp.rowCenter.txt          --plotHeight 20   --plotWidth 20      >> $output_sub1/A.7A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/A.1-results.10000bp.npz     --ntop 10000    --rowCenter                        -o $output_sub1/A.7B-top10000-plotPCA.10000bp.rowCenter.pdf       --outFileNameData  $output_sub1/A.7B-top10000-plotPCA.10000bp.rowCenter.txt         --plotHeight 20   --plotWidth 20      >> $output_sub1/A.7B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/A.1-results.10000bp.npz     --ntop 0        --rowCenter                        -o $output_sub1/A.7C-topAll-plotPCA.10000bp.rowCenter.pdf         --outFileNameData  $output_sub1/A.7C-topAll-plotPCA.10000bp.rowCenter.txt           --plotHeight 20   --plotWidth 20      >> $output_sub1/A.7C-plotPCA-runLog.txt   2>&1");                               

     $output_sub1 = "$output_g/multiBigwigSummary_bin/5000bp";   
     &myMakeDir($output_sub1);
    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 5000   --numberOfProcessors max/2    --verbose    --outRawCounts $output_sub1/B.1-RawCounts.5000bp-Bin.txt   --outFileName $output_sub1/B.1-results.5000bp.npz    >> $output_sub1/B.1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/B.1-results.5000bp.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/B.2-Correlation.heatmap.pearson.5000bp.pdf       --outFileCorMatrix $output_sub1/B.2-Correlation.heatmap.pearson.5000bp.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/B.2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/B.1-results.5000bp.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/B.3-Correlation.heatmap.spearman.5000bp.pdf      --outFileCorMatrix $output_sub1/B.3-Correlation.heatmap.spearman.5000bp.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/B.3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/B.1-results.5000bp.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/B.4-Correlation.scatterplot.pearson.5000bp.pdf   --outFileCorMatrix $output_sub1/B.4-Correlation.scatterplot.pearson.5000bp.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/B.4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/B.1-results.5000bp.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/B.5-Correlation.scatterplot.spearman.5000bp.pdf  --outFileCorMatrix $output_sub1/B.5-Correlation.scatterplot.spearman.5000bp.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/B.5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/B.1-results.5000bp.npz                                                        -o $output_sub1/B.6A-top1000-plotPCA.5000bp.pdf                  --outFileNameData  $output_sub1/B.6A-top1000-plotPCA.5000bp.txt                    --plotHeight 20   --plotWidth 20      >> $output_sub1/B.6A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/B.1-results.5000bp.npz     --ntop 10000                                       -o $output_sub1/B.6B-top10000-plotPCA.5000bp.pdf                 --outFileNameData  $output_sub1/B.6B-top10000-plotPCA.5000bp.txt                   --plotHeight 20   --plotWidth 20      >> $output_sub1/B.6B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/B.1-results.5000bp.npz     --ntop 0                                           -o $output_sub1/B.6C-topAll-plotPCA.5000bp.pdf                   --outFileNameData  $output_sub1/B.6C-topAll-plotPCA.5000bp.txt                     --plotHeight 20   --plotWidth 20      >> $output_sub1/B.6C-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/B.1-results.5000bp.npz                     --rowCenter                        -o $output_sub1/B.7A-top1000-plotPCA.5000bp.rowCenter.pdf        --outFileNameData  $output_sub1/B.7A-top1000-plotPCA.5000bp.rowCenter.txt          --plotHeight 20   --plotWidth 20      >> $output_sub1/B.7A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/B.1-results.5000bp.npz     --ntop 10000    --rowCenter                        -o $output_sub1/B.7B-top10000-plotPCA.5000bp.rowCenter.pdf       --outFileNameData  $output_sub1/B.7B-top10000-plotPCA.5000bp.rowCenter.txt         --plotHeight 20   --plotWidth 20      >> $output_sub1/B.7B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/B.1-results.5000bp.npz     --ntop 0        --rowCenter                        -o $output_sub1/B.7C-topAll-plotPCA.5000bp.rowCenter.pdf         --outFileNameData  $output_sub1/B.7C-topAll-plotPCA.5000bp.rowCenter.txt           --plotHeight 20   --plotWidth 20      >> $output_sub1/B.7C-plotPCA-runLog.txt   2>&1");                               

     $output_sub1 = "$output_g/multiBigwigSummary_bin/1000bp";   
     &myMakeDir($output_sub1);
    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 1000   --numberOfProcessors max/2    --verbose    --outRawCounts $output_sub1/C.1-RawCounts.1000bp-Bin.txt   --outFileName $output_sub1/C.1-results.1000bp.npz    >> $output_sub1/C.1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/C.1-results.1000bp.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/C.2-Correlation.heatmap.pearson.1000bp.pdf       --outFileCorMatrix $output_sub1/C.2-Correlation.heatmap.pearson.1000bp.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/C.2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/C.1-results.1000bp.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/C.3-Correlation.heatmap.spearman.1000bp.pdf      --outFileCorMatrix $output_sub1/C.3-Correlation.heatmap.spearman.1000bp.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/C.3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/C.1-results.1000bp.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/C.4-Correlation.scatterplot.pearson.1000bp.pdf   --outFileCorMatrix $output_sub1/C.4-Correlation.scatterplot.pearson.1000bp.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/C.4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/C.1-results.1000bp.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/C.5-Correlation.scatterplot.spearman.1000bp.pdf  --outFileCorMatrix $output_sub1/C.5-Correlation.scatterplot.spearman.1000bp.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/C.5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/C.1-results.1000bp.npz                                                        -o $output_sub1/C.6A-top1000-plotPCA.1000bp.pdf                  --outFileNameData  $output_sub1/C.6A-top1000-plotPCA.1000bp.txt                    --plotHeight 20   --plotWidth 20      >> $output_sub1/C.6A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/C.1-results.1000bp.npz     --ntop 10000                                       -o $output_sub1/C.6B-top10000-plotPCA.1000bp.pdf                 --outFileNameData  $output_sub1/C.6B-top10000-plotPCA.1000bp.txt                   --plotHeight 20   --plotWidth 20      >> $output_sub1/C.6B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/C.1-results.1000bp.npz     --ntop 0                                           -o $output_sub1/C.6C-topAll-plotPCA.1000bp.pdf                   --outFileNameData  $output_sub1/C.6C-topAll-plotPCA.1000bp.txt                     --plotHeight 20   --plotWidth 20      >> $output_sub1/C.6C-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/C.1-results.1000bp.npz                     --rowCenter                        -o $output_sub1/C.7A-top1000-plotPCA.1000bp.rowCenter.pdf        --outFileNameData  $output_sub1/C.7A-top1000-plotPCA.1000bp.rowCenter.txt          --plotHeight 20   --plotWidth 20      >> $output_sub1/C.7A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/C.1-results.1000bp.npz     --ntop 10000    --rowCenter                        -o $output_sub1/C.7B-top10000-plotPCA.1000bp.rowCenter.pdf       --outFileNameData  $output_sub1/C.7B-top10000-plotPCA.1000bp.rowCenter.txt         --plotHeight 20   --plotWidth 20      >> $output_sub1/C.7B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/C.1-results.1000bp.npz     --ntop 0        --rowCenter                        -o $output_sub1/C.7C-topAll-plotPCA.1000bp.rowCenter.pdf         --outFileNameData  $output_sub1/C.7C-topAll-plotPCA.1000bp.rowCenter.txt           --plotHeight 20   --plotWidth 20      >> $output_sub1/C.7C-plotPCA-runLog.txt   2>&1");                               

    $output_sub1 = "$output_g/multiBigwigSummary_bin/500bp";   
    &myMakeDir($output_sub1);
    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 500   --numberOfProcessors max/2    --verbose    --outRawCounts $output_sub1/D.1-RawCounts.500bp-Bin.txt   --outFileName $output_sub1/D.1-results.500bp.npz    >> $output_sub1/D.1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/D.1-results.500bp.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/D.2-Correlation.heatmap.pearson.500bp.pdf       --outFileCorMatrix $output_sub1/D.2-Correlation.heatmap.pearson.500bp.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/D.2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/D.1-results.500bp.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/D.3-Correlation.heatmap.spearman.500bp.pdf      --outFileCorMatrix $output_sub1/D.3-Correlation.heatmap.spearman.500bp.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/D.3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/D.1-results.500bp.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/D.4-Correlation.scatterplot.pearson.500bp.pdf   --outFileCorMatrix $output_sub1/D.4-Correlation.scatterplot.pearson.500bp.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/D.4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/D.1-results.500bp.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/D.5-Correlation.scatterplot.spearman.500bp.pdf  --outFileCorMatrix $output_sub1/D.5-Correlation.scatterplot.spearman.500bp.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/D.5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/D.1-results.500bp.npz                                                        -o $output_sub1/D.6A-top1000-plotPCA.500bp.pdf                  --outFileNameData  $output_sub1/D.6A-top1000-plotPCA.500bp.txt                    --plotHeight 20   --plotWidth 20      >> $output_sub1/D.6A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/D.1-results.500bp.npz     --ntop 10000                                       -o $output_sub1/D.6B-top10000-plotPCA.500bp.pdf                 --outFileNameData  $output_sub1/D.6B-top10000-plotPCA.500bp.txt                   --plotHeight 20   --plotWidth 20      >> $output_sub1/D.6B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/D.1-results.500bp.npz     --ntop 0                                           -o $output_sub1/D.6C-topAll-plotPCA.500bp.pdf                   --outFileNameData  $output_sub1/D.6C-topAll-plotPCA.500bp.txt                     --plotHeight 20   --plotWidth 20      >> $output_sub1/D.6C-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/D.1-results.500bp.npz                     --rowCenter                        -o $output_sub1/D.7A-top1000-plotPCA.500bp.rowCenter.pdf        --outFileNameData  $output_sub1/D.7A-top1000-plotPCA.500bp.rowCenter.txt          --plotHeight 20   --plotWidth 20      >> $output_sub1/D.7A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/D.1-results.500bp.npz     --ntop 10000    --rowCenter                        -o $output_sub1/D.7B-top10000-plotPCA.500bp.rowCenter.pdf       --outFileNameData  $output_sub1/D.7B-top10000-plotPCA.500bp.rowCenter.txt         --plotHeight 20   --plotWidth 20      >> $output_sub1/D.7B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/D.1-results.500bp.npz     --ntop 0        --rowCenter                        -o $output_sub1/D.7C-topAll-plotPCA.500bp.rowCenter.pdf         --outFileNameData  $output_sub1/D.7C-topAll-plotPCA.500bp.rowCenter.txt           --plotHeight 20   --plotWidth 20      >> $output_sub1/D.7C-plotPCA-runLog.txt   2>&1");                               

    $output_sub1 = "$output_g/multiBigwigSummary_bin/100bp";   
    &myMakeDir($output_sub1);
    system("multiBigwigSummary bins   --bwfiles @BigWigfiles_g2    --smartLabels    --binSize 100   --numberOfProcessors max/2    --verbose    --outRawCounts $output_sub1/E.1-RawCounts.100bp-Bin.txt   --outFileName $output_sub1/E.1-results.100bp.npz    >> $output_sub1/E.1-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/E.1-results.100bp.npz    --whatToPlot heatmap        --corMethod pearson     -o $output_sub1/E.2-Correlation.heatmap.pearson.100bp.pdf       --outFileCorMatrix $output_sub1/E.2-Correlation.heatmap.pearson.100bp.txt         --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/E.2-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/E.1-results.100bp.npz    --whatToPlot heatmap        --corMethod spearman    -o $output_sub1/E.3-Correlation.heatmap.spearman.100bp.pdf      --outFileCorMatrix $output_sub1/E.3-Correlation.heatmap.spearman.100bp.txt        --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/E.3-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/E.1-results.100bp.npz    --whatToPlot scatterplot    --corMethod pearson     -o $output_sub1/E.4-Correlation.scatterplot.pearson.100bp.pdf   --outFileCorMatrix $output_sub1/E.4-Correlation.scatterplot.pearson.100bp.txt     --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/E.4-runLog.txt   2>&1");                               
    system("plotCorrelation -in $output_sub1/E.1-results.100bp.npz    --whatToPlot scatterplot    --corMethod spearman    -o $output_sub1/E.5-Correlation.scatterplot.spearman.100bp.pdf  --outFileCorMatrix $output_sub1/E.5-Correlation.scatterplot.spearman.100bp.txt    --plotHeight 20   --plotWidth 20   --plotNumbers   --skipZeros  --removeOutliers     >> $output_sub1/E.5-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/E.1-results.100bp.npz                                                        -o $output_sub1/E.6A-top1000-plotPCA.100bp.pdf                  --outFileNameData  $output_sub1/E.6A-top1000-plotPCA.100bp.txt                    --plotHeight 20   --plotWidth 20      >> $output_sub1/E.6A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/E.1-results.100bp.npz     --ntop 10000                                       -o $output_sub1/E.6B-top10000-plotPCA.100bp.pdf                 --outFileNameData  $output_sub1/E.6B-top10000-plotPCA.100bp.txt                   --plotHeight 20   --plotWidth 20      >> $output_sub1/E.6B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/E.1-results.100bp.npz     --ntop 0                                           -o $output_sub1/E.6C-topAll-plotPCA.100bp.pdf                   --outFileNameData  $output_sub1/E.6C-topAll-plotPCA.100bp.txt                     --plotHeight 20   --plotWidth 20      >> $output_sub1/E.6C-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/E.1-results.100bp.npz                     --rowCenter                        -o $output_sub1/E.7A-top1000-plotPCA.100bp.rowCenter.pdf        --outFileNameData  $output_sub1/E.7A-top1000-plotPCA.100bp.rowCenter.txt          --plotHeight 20   --plotWidth 20      >> $output_sub1/E.7A-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/E.1-results.100bp.npz     --ntop 10000    --rowCenter                        -o $output_sub1/E.7B-top10000-plotPCA.100bp.rowCenter.pdf       --outFileNameData  $output_sub1/E.7B-top10000-plotPCA.100bp.rowCenter.txt         --plotHeight 20   --plotWidth 20      >> $output_sub1/E.7B-plotPCA-runLog.txt   2>&1");                               
    system("plotPCA         -in $output_sub1/E.1-results.100bp.npz     --ntop 0        --rowCenter                        -o $output_sub1/E.7C-topAll-plotPCA.100bp.rowCenter.pdf         --outFileNameData  $output_sub1/E.7C-topAll-plotPCA.100bp.rowCenter.txt           --plotHeight 20   --plotWidth 20      >> $output_sub1/E.7C-plotPCA-runLog.txt   2>&1");                               

}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END







