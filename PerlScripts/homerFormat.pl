#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;

## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "1-renamed"   
my $output_g = '';  ## such as "2-formatted"

my $HELP = "  perl  homerFormat.pl -in 1-renamed   -out 2-formatted  ";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '1-renamed';     ## This is only an initialization value or suggesting value, not default value.
$output_g = '2-formatted';   ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "    -in   -out  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  homerFormat.pl  -help' \n";
    exit 0;
}

## Get Arguments 
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
 
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";


sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}

&myMakeDir("$output_g");

opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
###################################################################################################################################################################################################








###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] =~ m/\.txt$/;
        my $temp = $inputFiles_g[$i];
        say   "\t...... $temp " ;
        open(INPUT1_FH,   "<",   "$input_g/$temp"  )      or   die "$!"; 
        open(OUTPUT1_FH,  ">",   "$output_g/$temp" )      or   die "$!"; 
        my @lines1 = <INPUT1_FH>; 
        print  OUTPUT1_FH   "chr\tstart\tend\tstrand\tname\tstrength\n";

        for (my $i=1; $i<=$#lines1; $i++) {
            my $temp1 = $lines1[$i];
            $temp1 =~ m/^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s*/ or die;
            my $chr   = $2;
            my $start = $3;
            my $end   = $4 + 1;
            my $strand = "*";
            my $name   = $temp."_".$i;
            my $strength = "1";
            print  OUTPUT1_FH   "$chr\t$start\t$end\t$strand\t$name\t$strength\n";
        }

}


}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
