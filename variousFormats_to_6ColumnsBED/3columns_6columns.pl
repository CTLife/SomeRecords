#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;
## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
## Use the 1st to 3rd columns to genarate 6 columns BED files.
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as dir "10-pooled"  
my $outDir_g = '';  ## such as dir "10-bed6cols"    
my $name_g   = '';  ## Name of genomic regions
my $HELP     = "  perl  3columns_6columns.pl  -in xxx  -out xxx  -name xxx ";
## perl  3columns_6columns.pl   -in 10-pooled     -out 10-bed6cols   -name DMRs

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '10-pooled';     ## This is only an initialization value or suggesting value, not default value.
$outDir_g = '10-bed6cols';   ## This is only an initialization value or suggesting value, not default value.
$name_g   = 'DMRs';          ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "    -in   -out  -name  ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl  3columns_6columns.pl  -help' \n";
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
&myMakeDir($outDir_g);
&myMakeDir("$outDir_g/BED");
&myMakeDir("$outDir_g/BED_withHeader");

opendir(my $DH_input_g, $input_g)  ||  die;
my @inputFiles_g = readdir($DH_input_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
{
say   "\n\n\n\n\n\n##################################################################################################";

for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        next unless $inputFiles_g[$i] !~ m/^[.]/;
        next unless $inputFiles_g[$i] !~ m/[~]$/;
        next unless $inputFiles_g[$i] =~ m/\.bed$/;
        my $temp = $inputFiles_g[$i];
        say   "\t...... $temp " ;
        $temp =~ s/.bed$// or die;
        open(INPUT1_FH,   "<",   "$input_g/$temp.bed" )      or   die "$!"; 
        open(OUTPUT1_FH,  ">",   "$outDir_g/BED/$temp.bed" )      or   die "$!"; 
        open(OUTPUT2_FH,  ">",   "$outDir_g/BED_withHeader/$temp.bed" )      or   die "$!"; 

        my @lines1 = <INPUT1_FH>; 
        print  OUTPUT2_FH   "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n";

        for (my $j=0; $j<=$#lines1; $j++) {
            my $temp1 = $lines1[$j];
            $temp1 =~ m/^(\S+)\s+(\d+)\s+(\d+)\s*/ or die  "$! \n\n $temp1 \n\n ";
            my $chr   = $1;
            my $start = $2;
            my $end   = $3;
            my $strand = "*";
            my $name   = $name_g."_".$temp."_".$j;
            my $strength = "1";
            print  OUTPUT1_FH   "$chr\t$start\t$end\t$name\t$strength\t$strand\n";
            print  OUTPUT2_FH   "$chr\t$start\t$end\t$name\t$strength\t$strand\n";

        }

}


}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
