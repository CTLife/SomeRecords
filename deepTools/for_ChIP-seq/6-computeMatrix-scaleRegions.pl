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
my $output_g = '';  ## such as "7-computeMatrix-scaleRegions/H3K4me3"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Usage:
               perl  6-computeMatrix-scaleRegions.pl    [-version]    [-help]      [-in inputDir]     [-bed bedDir]     [-out outDir]
        For instance:
               perl  6-computeMatrix-scaleRegions.pl    -in 3-subtractInput/H3K4me3   -bed GenomicRegions/3_Enhancers/allPups    -out 7-computeMatrix-scaleRegions/H3K4me3    > 6-computeMatrix-scaleRegions.H3K4me3.runLog

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
$input_g    = '3-subtractInput/H3K4me3';                ## This is only an initialization value or suggesting value, not default value.
$BEDdir_g   = '3_Enhancers/allPups';                    ## This is only an initialization value or suggesting value, not default value.
$output_g   = '7-computeMatrix-scaleRegions/H3K4me3';   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help   -genome   -bed   -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  6-computeMatrix-scaleRegions.pl  -help' \n";
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
say   "Using computeMatrix  reference-point, plotProfile and plotHeatmap ......";

for ( my $i=0; $i<=$#BED_Files_g; $i++ ) { 
    next unless $BED_Files_g[$i] =~ m/\.bed$/;
    next unless $BED_Files_g[$i] !~ m/^[.]/;
    next unless $BED_Files_g[$i] !~ m/[~]$/;
    for ( my $j=0; $j<=$#BigWigfiles_g; $j++ ) {
         my $output_sub1 = "$output_g/$BED_Files_g[$i]/$BigWigfiles_g[$j]";   
         &myMakeDir($output_sub1);
         system("computeMatrix    scale-regions    --scoreFileName $input_g/$BigWigfiles_g[$j]    --regionsFileName $BEDdir_g/$BED_Files_g[$i]     --smartLabels   --numberOfProcessors max/2    --verbose    -o $output_sub1/1.gzipped-matrix.gz    --outFileNameMatrix $output_sub1/1.matrix.txt   --outFileSortedRegions  $output_sub1/1.Regions.txt    --regionBodyLength  1000   --upstream 1000   --downstream 1000  --binSize 50    --missingDataAsZero     >> $output_sub1/1.runLog.txt   2>&1");                               

         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2A-figure.plotProfile.lines.C=1.pdf   --outFileSortedRegions $output_sub1/2A-regions.plotProfile.lines.C=1.txt   --outFileNameData  $output_sub1/2A-Profile.plotProfile.lines.C=1.tab.txt                 --averageType mean  --plotType lines     --verbose     >> $output_sub1/2A.plotProfile.lines.C=1.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2B-figure.plotProfile.lines.C=2.pdf   --outFileSortedRegions $output_sub1/2B-regions.plotProfile.lines.C=2.txt   --outFileNameData  $output_sub1/2B-Profile.plotProfile.lines.C=2.tab.txt   --kmeans 2    --averageType mean  --plotType lines     --verbose     >> $output_sub1/2B.plotProfile.lines.C=2.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2C-figure.plotProfile.lines.C=3.pdf   --outFileSortedRegions $output_sub1/2C-regions.plotProfile.lines.C=3.txt   --outFileNameData  $output_sub1/2C-Profile.plotProfile.lines.C=3.tab.txt   --kmeans 3    --averageType mean  --plotType lines     --verbose     >> $output_sub1/2C.plotProfile.lines.C=3.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2D-figure.plotProfile.lines.C=4.pdf   --outFileSortedRegions $output_sub1/2D-regions.plotProfile.lines.C=4.txt   --outFileNameData  $output_sub1/2D-Profile.plotProfile.lines.C=4.tab.txt   --kmeans 4    --averageType mean  --plotType lines     --verbose     >> $output_sub1/2D.plotProfile.lines.C=4.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2E-figure.plotProfile.lines.C=5.pdf   --outFileSortedRegions $output_sub1/2E-regions.plotProfile.lines.C=5.txt   --outFileNameData  $output_sub1/2E-Profile.plotProfile.lines.C=5.tab.txt   --kmeans 5    --averageType mean  --plotType lines     --verbose     >> $output_sub1/2E.plotProfile.lines.C=5.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2F-figure.plotProfile.lines.C=7.pdf   --outFileSortedRegions $output_sub1/2F-regions.plotProfile.lines.C=7.txt   --outFileNameData  $output_sub1/2F-Profile.plotProfile.lines.C=7.tab.txt   --kmeans 7    --averageType mean  --plotType lines     --verbose     >> $output_sub1/2F.plotProfile.lines.C=7.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/2G-figure.plotProfile.lines.C=9.pdf   --outFileSortedRegions $output_sub1/2G-regions.plotProfile.lines.C=9.txt   --outFileNameData  $output_sub1/2G-Profile.plotProfile.lines.C=9.tab.txt   --kmeans 9    --averageType mean  --plotType lines     --verbose     >> $output_sub1/2G.plotProfile.lines.C=9.runLog.txt   2>&1 "  );

         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3A-figure.plotProfile.se.C=1.pdf   --outFileSortedRegions $output_sub1/3A-regions.plotProfile.se.C=1.txt   --outFileNameData  $output_sub1/3A-Profile.plotProfile.se.C=1.tab.txt                 --averageType mean  --plotType se     --verbose     >> $output_sub1/3A.plotProfile.se.C=1.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3B-figure.plotProfile.se.C=2.pdf   --outFileSortedRegions $output_sub1/3B-regions.plotProfile.se.C=2.txt   --outFileNameData  $output_sub1/3B-Profile.plotProfile.se.C=2.tab.txt   --kmeans 2    --averageType mean  --plotType se     --verbose     >> $output_sub1/3B.plotProfile.se.C=2.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3C-figure.plotProfile.se.C=3.pdf   --outFileSortedRegions $output_sub1/3C-regions.plotProfile.se.C=3.txt   --outFileNameData  $output_sub1/3C-Profile.plotProfile.se.C=3.tab.txt   --kmeans 3    --averageType mean  --plotType se     --verbose     >> $output_sub1/3C.plotProfile.se.C=3.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3D-figure.plotProfile.se.C=4.pdf   --outFileSortedRegions $output_sub1/3D-regions.plotProfile.se.C=4.txt   --outFileNameData  $output_sub1/3D-Profile.plotProfile.se.C=4.tab.txt   --kmeans 4    --averageType mean  --plotType se     --verbose     >> $output_sub1/3D.plotProfile.se.C=4.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3E-figure.plotProfile.se.C=5.pdf   --outFileSortedRegions $output_sub1/3E-regions.plotProfile.se.C=5.txt   --outFileNameData  $output_sub1/3E-Profile.plotProfile.se.C=5.tab.txt   --kmeans 5    --averageType mean  --plotType se     --verbose     >> $output_sub1/3E.plotProfile.se.C=5.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3F-figure.plotProfile.se.C=7.pdf   --outFileSortedRegions $output_sub1/3F-regions.plotProfile.se.C=7.txt   --outFileNameData  $output_sub1/3F-Profile.plotProfile.se.C=7.tab.txt   --kmeans 7    --averageType mean  --plotType se     --verbose     >> $output_sub1/3F.plotProfile.se.C=7.runLog.txt   2>&1 "  );
         system("plotProfile   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/3G-figure.plotProfile.se.C=9.pdf   --outFileSortedRegions $output_sub1/3G-regions.plotProfile.se.C=9.txt   --outFileNameData  $output_sub1/3G-Profile.plotProfile.se.C=9.tab.txt   --kmeans 9    --averageType mean  --plotType se     --verbose     >> $output_sub1/3G.plotProfile.se.C=9.runLog.txt   2>&1 "  );

         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5A-figure.plotHeatmap.C=1.pdf   --outFileSortedRegions $output_sub1/5A-regions.plotHeatmap.C=1.txt   --outFileNameMatrix  $output_sub1/5A-Matrix.plotHeatmap.C=1.tab.txt                 --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5A.plotHeatmap.C=1.runLog.txt   2>&1 "  );
         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5B-figure.plotHeatmap.C=2.pdf   --outFileSortedRegions $output_sub1/5B-regions.plotHeatmap.C=2.txt   --outFileNameMatrix  $output_sub1/5B-Matrix.plotHeatmap.C=2.tab.txt   --kmeans 2    --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5B.plotHeatmap.C=2.runLog.txt   2>&1 "  );
         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5C-figure.plotHeatmap.C=3.pdf   --outFileSortedRegions $output_sub1/5C-regions.plotHeatmap.C=3.txt   --outFileNameMatrix  $output_sub1/5C-Matrix.plotHeatmap.C=3.tab.txt   --kmeans 3    --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5C.plotHeatmap.C=3.runLog.txt   2>&1 "  );
         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5D-figure.plotHeatmap.C=4.pdf   --outFileSortedRegions $output_sub1/5D-regions.plotHeatmap.C=4.txt   --outFileNameMatrix  $output_sub1/5D-Matrix.plotHeatmap.C=4.tab.txt   --kmeans 4    --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5D.plotHeatmap.C=4.runLog.txt   2>&1 "  );
         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5E-figure.plotHeatmap.C=5.pdf   --outFileSortedRegions $output_sub1/5E-regions.plotHeatmap.C=5.txt   --outFileNameMatrix  $output_sub1/5E-Matrix.plotHeatmap.C=5.tab.txt   --kmeans 5    --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5E.plotHeatmap.C=5.runLog.txt   2>&1 "  );
         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5F-figure.plotHeatmap.C=7.pdf   --outFileSortedRegions $output_sub1/5F-regions.plotHeatmap.C=7.txt   --outFileNameMatrix  $output_sub1/5F-Matrix.plotHeatmap.C=7.tab.txt   --kmeans 7    --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5F.plotHeatmap.C=7.runLog.txt   2>&1 "  );
         system("plotHeatmap   --matrixFile $output_sub1/1.gzipped-matrix.gz   --outFileName $output_sub1/5G-figure.plotHeatmap.C=9.pdf   --outFileSortedRegions $output_sub1/5G-regions.plotHeatmap.C=9.txt   --outFileNameMatrix  $output_sub1/5G-Matrix.plotHeatmap.C=9.tab.txt   --kmeans 9    --sortRegions descend    --sortUsing mean     --verbose     >> $output_sub1/5G.plotHeatmap.C=9.runLog.txt   2>&1 "  );

    }
}
###################################################################################################################################################################################################    





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";


 

  
## END







