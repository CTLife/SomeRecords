#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1-raw/mergedAll.snp", global variable.
my $output_g = '';  ## such as "2-sepAGCT/mergedAll.snp", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Separate each line based on the 3rd column "ref" in the varscan2 output files.

        Usage:
               perl  2-sepAGCT_varscan.pl    [-version]    [-help]     [-in inputFile]    [-out outFile]
        For instance:
               perl  2-sepAGCT_varscan.pl    -in 1-raw/mergedAll.snp   -out 2-sepAGCT/mergedAll.snp    > 2-sepAGCT_varscan.runLog.txt

        Optional arguments:
        -version        Show version number of this program and exit.
        -help           Show this help message and exit.

        Required arguments:
        -in inputFile        "inputFile" is the name of your input file which is the output of VarScan.  (no default)
        -out outFile         "outFile" is the name of the running results of this script.  (no default)
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    Version 0.2,  2020-08-01.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '1-raw/mergedAll.snp';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '2-sepAGCT/mergedAll.snp';     ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  2-sepAGCT_varscan.pl  -help' \n";
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
                Input   File:  $input_g
                Output  Path:  $output_g
        ##########################################################
\n";
}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}
&myMakeDir($output_g);

open(fileA, ">", "$output_g/A.txt")  or  die;
open(fileG, ">", "$output_g/G.txt")  or  die;
open(fileC, ">", "$output_g/C.txt")  or  die;
open(fileT, ">", "$output_g/T.txt")  or  die;
open(fileN, ">", "$output_g/N.txt")  or  die;

open(normal_bed_fileA, ">", "$output_g/Znormal_A.bed")  or  die;
open(normal_bed_fileG, ">", "$output_g/Znormal_G.bed")  or  die;
open(normal_bed_fileC, ">", "$output_g/Znormal_C.bed")  or  die;
open(normal_bed_fileT, ">", "$output_g/Znormal_T.bed")  or  die;
open(normal_bed_fileN, ">", "$output_g/Znormal_N.bed")  or  die;

open(Selected_bed_fileA, ">", "$output_g/Selected.A.bed")  or  die;
open(Selected_bed_fileG, ">", "$output_g/Selected.G.bed")  or  die;
open(Selected_bed_fileC, ">", "$output_g/Selected.C.bed")  or  die;
open(Selected_bed_fileT, ">", "$output_g/Selected.T.bed")  or  die;

my $minReads = 3;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
open(my $FH_input_g, $input_g)  ||  die;
my $line1 = <$FH_input_g>;
print fileA $line1;
print fileG $line1;
print fileC $line1;
print fileT $line1;
print fileN $line1;

my $n  = 0;
my $n1 = 0;
my $n2 = 0;
my $n3 = 0;
my $n4 = 0;
my $n5 = 0;

while(my $line=<$FH_input_g>) {
    $line =~ m/^(\S+)\t(\d+)\t([A-Z])\t(\S+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\S+)\n$/ or die "\n\n##$line##\n\n\n";
    my $Chrom    = $1;	
    my $Position = $2;	
    my $Ref      = $3;	
    my $Cons     = $4;	
    my $Reads1   = $5;	
    my $Reads2   = $6;	
    my $VarFreq  = $7;	
    my $Strands1 = $8;	
    my $Strands2 = $9;	
    my $Qual1    = $10;	
    my $Qual2    = $11;	
    my $Pvalue   = $12;	
    my $MapQual1 = $13;	
    my $MapQual2 = $14;	
    my $Reads1Plus   = $15;	
    my $Reads1Minus  = $16;	
    my $Reads2Plus   = $17;	
    my $Reads2Minus  = $18;	
    my $VarAllele    = $19;	

    my $start  = $Position - 1; 
    my $end    = $Position; 
    my $name1  = "$Chrom...$Position...$Ref...$Cons...$Reads1...$Reads2...$VarFreq...$Strands1...$Strands2...$Qual1...$Qual2...$Pvalue";
    my $name2  = "$MapQual1...$MapQual2...$Reads1Plus...$Reads1Minus...$Reads2Plus...$Reads2Minus...$VarAllele";                               
    my $name   = "$name1...$name2"; 
    $VarFreq =~ s/%// or die;

    $n++;
    given($Ref) {
        	when("A") { print fileA $line; 
                            print normal_bed_fileA "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; 
                            if( ($VarFreq >= 0) and ($Reads1+$Reads2 >= $minReads) and ($VarAllele eq "G") ) {print Selected_bed_fileA "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; }                                
                            $n1++;
                          }
		when("G") { print fileG $line; 
                            print normal_bed_fileG "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; 
                            if( ($VarFreq >= 0) and ($Reads1+$Reads2 >= $minReads) and ($VarAllele eq "A") ) {print Selected_bed_fileG "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; }                                 
                            $n2++;
                          }
		when("C") { print fileC $line; 
                            print normal_bed_fileC "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; 
                            if( ($VarFreq >= 0) and ($Reads1+$Reads2 >= $minReads) and ($VarAllele eq "T") ) {print Selected_bed_fileC "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; }                                  
                            $n3++;
                          }
		when("T") { print fileT $line; 
                            print normal_bed_fileT "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n";  
                            if( ($VarFreq >= 0) and ($Reads1+$Reads2 >= $minReads) and ($VarAllele eq "C") ) {print Selected_bed_fileT "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n"; }                                
                            $n4++;
                          }
		default   { print fileN $line; 
                            print normal_bed_fileN "$Chrom\t$start\t$end\t$name\t$VarFreq\t.\t$Ref\t$VarAllele\t$Reads1\t$Reads2\n";  
                            $n5++;
                          }
        }
    
}



my $m = $n1+$n2+$n3+$n4+$n5;
if($n == $m){
         print "OK:\n";
         print "$n == $m; $n1,$n2,$n3,$n4,$n5\n\n";
}else{
         print "Wrong:\n";
         print "$n != $m; $n1,$n2,$n3,$n4,$n5\n\n";     
}


#####
print  "Done!!!\n\n\n\n\n";













