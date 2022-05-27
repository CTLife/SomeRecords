#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "8-select/GY1", global variable.
my $output_g = '';  ## such as "10-mergeNearest/GY1", global variable.
{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------

        Usage:
               perl  10-mergeNearest.pl    [-version]    [-help]     [-in inputFile]    [-out outFile]
        For instance:
               perl  10-mergeNearest.pl    -in 8-select/GY1   -out 10-mergeNearest/GY1    > 10-mergeNearest.GY1.A.runLog.txt

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
$input_g  = '8-select/GY1';         ## This is only an initialization value or suggesting value, not default value.
$output_g = '10-mergeNearest/GY1';        ## This is only an initialization value or suggesting value, not default value.

## Available Arguments
my $available = "   -version    -help    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  10-mergeNearest.pl  -help' \n";
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
&myMakeDir( "$output_g.selected" );
###################################################################################################################################################################################################





###################################################################################################################################################################################################
opendir(my $FH_input_g, $input_g)  ||  die;
my @files = readdir($FH_input_g);

for(my $i=0; $i<=$#files; $i++) {
    next unless $files[$i] =~ m/\.bed$/; 
    my $input1 = "$input_g/$files[$i]";
    system(  "bedtools sort  -i $input1  > $input1.sorted" );
    
    my $output1 = "$output_g/$files[$i]";
    my $output2 = "$output_g.selected/$files[$i]";
    open(FH_input,  "<",  "$input1.sorted")   or  die;
    open(FH_output1, ">",   $output1)  or  die;
    open(FH_output2, ">",   $output2)  or  die;
    my @lines1 = <FH_input>;
    
    my $merged_chr     = "";
    my $merged_start   = "";
    my $merged_end     = "";
    my $merged_name    = "";
    my $merged_num     = 0;
    my $merged_strand  = ".";
    
    for(my $j=0; $j<=$#lines1; $j++) {
          $lines1[$j] =~ m/^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+/ or die  "## $lines1[$j] ##\n"; 
          my $chr   = $1;
          my $start = $2;
          my $end   = $3;
          my $name  = $4;
          
          if($j == 0 ) { 
    		 $merged_chr     = $chr;
    		 $merged_start   = $start;
   		 $merged_end     = $end;                           
          }
          
          my $diff = $end - $merged_end ; 
          if(  ($merged_chr eq $chr) and ( $diff <= 20 )  ) {  
          	$merged_end  = $end ;  
          	$merged_name = $merged_name."......".$name;   
          	$merged_num  = $merged_num+1;
          }else{
                 $merged_end = $merged_end + 1; 
                 print FH_output1 "$merged_chr\t$merged_start\t$merged_end\t$merged_name\t$merged_num\t$merged_strand\n";
                 if($merged_num >= 5) {print FH_output2  "$merged_chr\t$merged_start\t$merged_end\t$merged_name\t$merged_num\t$merged_strand\n";}
                 ($merged_num >= 2) or die "\n## $merged_num ##\n\n";
                 $merged_chr     = $chr;
   		 $merged_start   = $start;
   		 $merged_end     = $end;
   		 $merged_name    = $name;
    		 $merged_num     = 1;
          }       
    }
        
}




#####
print  "Done!!!\n\n\n\n\n";













