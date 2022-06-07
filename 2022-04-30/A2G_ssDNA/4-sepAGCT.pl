#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "3-removeSNP/GY1/A.bed", global variable.
my $output_g = '';  ## such as "4-sepVar/GY1/A", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        Separate each line based on the variation bases in the varscan2 output files.

        Usage:
               perl  4-sepAGCT.pl    [-version]    [-help]     [-in inputFile]    [-out outFile]
        For instance:
               perl  4-sepAGCT.pl    -in 3-removeSNP/GY1/A.bed   -out 4-sepVar/GY1/A    > 4-sepAGCT.GY1.A.runLog.txt

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
$input_g  = '3-removeSNP/GY1/A.bed';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '4-sepVar/GY1/A';                ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  4-sepAGCT.pl  -help' \n";
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

open(fileA, ">", "$output_g/A.bed")  or  die;
open(fileG, ">", "$output_g/G.bed")  or  die;
open(fileC, ">", "$output_g/C.bed")  or  die;
open(fileT, ">", "$output_g/T.bed")  or  die;
open(fileN, ">", "$output_g/N.bed")  or  die;

my $minReads = 3;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
open(my $FH_input_g, $input_g)  ||  die;
my $line1 = <$FH_input_g>;
my $n  = 0;
my $n1 = 0;
my $n2 = 0;
my $n3 = 0;
my $n4 = 0;
my $n5 = 0;

while(my $line=<$FH_input_g>) {
    $line =~ m/\t(\S+)\t(\.)\t(\S)\t(\S+)\t(\d+)\t(\d+)\n$/ or die "\n\n##$line##\n\n\n";
    # my $VarFreq  = $1;    	
    # my $Ref      = $3;	
    my $VarAllele    = $4;    
    # my $Reads1   = $5;	
    # my $Reads2   = $6;	
    	
    $n++;
    given($VarAllele) {
        	when("A") { print fileA $line; 
                            $n1++;
                          }
		when("G") { print fileG $line; 
                            $n2++;
                          }
		when("C") { print fileC $line; 
                            $n3++;
                          }
		when("T") { print fileT $line; 
                            $n4++;
                          }
		default   { print fileN $line; 
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













