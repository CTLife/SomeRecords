#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.22;

## Perl5 version >= 5.22
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $genome_g = 'hg38';  ## such as "mm10", "ce11", "hg38".
my $input_g  = '';  ## such as "1-Coverage-CpG"
my $output_g = '';  ## such as "200-bigwig-5reads"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Welcome to use BSDA (BS-Seq Data Analyzer), version 0.9.0, 2017-10-01.
        BSDA is a Pipeline for Single-end and Paired-end BS-Seq Data Analysis by Integrating Lots of Softwares.                                                            

        Usage:
               perl  Cov2BW-5reads.pl    [-version]    [-help]     [-in inputDir]    [-out outDir]
        For instance:
               perl  Cov2BW-5reads.pl   -in 1-Coverage-CpG   -out 200-bigwig-5reads    > Cov2BW-5reads.runLog

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of input path that contains your Coverage files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        For more details about this pipeline and other NGS data analysis piplines, please visit https://github.com/CTLife/2ndGS_Pipelines

        Yong Peng @ Jie Qiao Lab, yongp@outlook.com, Key Laboratory of Assisted Reproduction at Third Hospital,
        Academy for Advanced Interdisciplinary Studies, and Peking-Tsinghua Center for Life Sciences (CLS), Peking University, China.
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    The Seventh Step of BSDA (BS-Seq Data Analyzer), version 0.9.0, 2017-10-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '1-Coverage-CpG';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '200-bigwig-5reads';      ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  Cov2BW-5reads.pl  -help' \n";
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
my $pattern_g    = "[-.0-9A-Za-z]+";
opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);

my $numCores_g   = 8;
###################################################################################################################################################################################################






###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Detecting Coverage files in input folder ......";
my @CovFiles_g = ();
{

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {     
    next unless $inputFiles_g[$i] =~ m/\.bismark.cov$/;
    next unless $inputFiles_g[$i] !~ m/^[.]/;
    next unless $inputFiles_g[$i] !~ m/[~]$/;
    say    "\t......$inputFiles_g[$i]"; 
    $inputFiles_g[$i] =~ m/^(\d+)_($pattern_g)_(Rep[1-9])\.bismark.cov$/  or  die;  
    $CovFiles_g[$#CovFiles_g+1] =  $inputFiles_g[$i];
    say        "\t\t\t\tBAM file:  $inputFiles_g[$i]\n";
}

say        "\t\t\t\tAll Coverage files:@CovFiles_g\n\n";
my $num1 = $#CovFiles_g + 1;
say         "\t\t\t\tThere are $num1 Coverage files.\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Convert coverrage file to bedGraph and BW ......";
&myMakeDir("$output_g/bedGraph"); 

for (my $i=0; $i<=$#CovFiles_g; $i++) {
    my $temp = $CovFiles_g[$i]; 
    $temp =~ s/\.bismark.cov$//  ||  die; 
    say   "\t......$CovFiles_g[$i]";
    open(INPUT_FH,   "<",   "$input_g/$temp.bismark.cov"  )             or die "$!"; 
    open(OUTPUT1_FH, ">",   "$output_g/bedGraph/$temp.bedGraph"  )      or die "$!"; 
    print  OUTPUT1_FH  "track type=bedGraph\n";
    my @lines1 = <INPUT_FH>; 

    for (my $i=0; $i<=$#lines1; $i++) {
        my $temp1 = $lines1[$i];
        $temp1 =~ m/^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s*/ or die;
        my $chrom1 = $1;  
        my $start1 = $2;  
        my $end1   = $3;  
        my $level1 = $4;  
        my $numC1  = $5;  
        my $numT1  = $6;  
        $chrom1 =~ m/^chr[\dXY]+$/     or  die "$temp1\n";
        ($level1>=0 and $level1<=100)  or  die "$temp1\n";
        my $total_num = $numC1 + $numT1;
        if($total_num >= 5) {print   OUTPUT1_FH  "$chrom1\t$start1\t$end1\t$level1\n";}
    }
    system("cat   $output_g/bedGraph/$temp.bedGraph   | sed '1d' | sort  -k1,1  -k2,2n     >  $output_g/bedGraph/$temp.sorted.bedGraph");  
    system("bedGraphToBigWig  $output_g/bedGraph/$temp.sorted.bedGraph   0-Others/Shortcuts/$genome_g/$genome_g.chrom.sizes      $output_g/$temp.bw ");          
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
