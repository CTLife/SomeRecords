#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;
## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "5_DMR_Annotation/2-Merged"  
my $outDir_g = '';  ## such as "5_DMR_Annotation/3-Formatted"    
my $name_g   = '';  ## Name of genomic regions
my $HELP     = "  perl  format.pl  -in xxx  -out xxx  -name xxx ";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '5_DMR_Annotation/2-Merged';     ## This is only an initialization value or suggesting value, not default value.
$outDir_g = '5_DMR_Annotation/3-Formatted';  ## This is only an initialization value or suggesting value, not default value.
$name_g   = '1_NC1_vs_NC2_Tiles';            ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "    -in   -out  -name  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  format.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-in'  })  { $input_g  = $args{'-in'    }; }else{say   "\n -in     is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out' })  { $outDir_g = $args{'-out'   }; }else{say   "\n -out    is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-name'})  { $name_g   = $args{'-name'  }; }else{say   "\n -name   is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$outDir_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$name_g   =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";


## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input           Path:  $input_g
                Out             Path:  $outDir_g
                Genomic Regions Name:  $name_g
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

my $myOut_1 = "$outDir_g/forChIPseeker";
my $myOut_2 = "$outDir_g/BED";
&myMakeDir("$myOut_1");
&myMakeDir("$myOut_2");

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
        $temp =~ s/.txt$// or die;
        open(INPUT1_FH,   "<",   "$input_g/$temp.txt" )      or   die "$!"; 
        open(OUTPUT1_FH,  ">",   "$myOut_1/$temp.txt" )      or   die "$!"; 
        open(OUTPUT2_FH,  ">",   "$myOut_2/$temp.bed" )      or   die "$!"; 

        my @lines1 = <INPUT1_FH>; 
        print  OUTPUT1_FH   "chr\tstart\tend\tstrand\tname\tstrength\n";

        for (my $j=4; $j<=$#lines1; $j++) {
            my $temp1 = $lines1[$j];
            $temp1 =~ m/^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s*/ or die  "$! \n\n $temp1 \n\n ";
            my $chr   = $2;
            my $start = $3-1;
            my $end   = $4;
            my $strand = "*";
            my $name   = $name_g."_".$temp."_".$j;
            my $strength = "1";

            print  OUTPUT1_FH   "$chr\t$start\t$end\t$strand\t$name\t$strength\n";
            print  OUTPUT2_FH   "$chr\t$start\t$end\t$name\n";

        }

}


}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
