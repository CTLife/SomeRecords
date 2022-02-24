#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my $input_g  = '';  ## such as "3_quan"   
my $output_g = '';  ## such as "4b_minusInput"

{
## Help Infromation
my $HELP = '
        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------     
        Usage:
               perl 4b_minusInput.pl    [-version]    [-help]    [-in inputDir]    [-out outDir]      
        For instance:
               nohup time  perl 4b_minusInput.pl   -in 3_quan   -out 4b_minusInput     > 4b_minusInput.runLog.txt  2>&1    &

        ----------------------------------------------------------------------------------------------------------
        Optional arguments:
        -version        Show version number of this program and exit.

        -help           Show this help message and exit.

        Required arguments:
        -in inputDir        "inputDir" is the name of input path that contains your files.  (no default)

        -out outDir         "outDir" is the name of output path that contains your running results of this step.  (no default)
        -----------------------------------------------------------------------------------------------------------

        ------------------------------------------------------------------------------------------------------------------------------------------------------
        ------------------------------------------------------------------------------------------------------------------------------------------------------
';

## Version Infromation
my $version = "    2022-01-27.";

## Keys and Values
if ($#ARGV   == -1)   { say  "\n$HELP\n";  exit 0;  }       ## when there are no any command argumants.
if ($#ARGV%2 ==  0)   { @ARGV = (@ARGV, "-help") ;  }       ## when the number of command argumants is odd.
my %args = @ARGV;

## Initialize  Variables
$input_g  = '3_quan';           ## This is only an initialization value or suggesting value, not default value.
$output_g = '4b_minusInput';    ## This is only an initialization value or suggesting value, not default value.

## Available Arguments 
my $available = "   -version    -help   -in   -out    ";
my $boole = 0;
while( my ($key, $value) = each %args ) {
    if ( ($key =~ m/^\-/) and ($available !~ m/\s$key\s/) ) {say    "\n\tCann't recognize $key";  $boole = 1; }
}
if($boole == 1) {
    say  "\tThe Command Line Arguments are wrong!";
    say  "\tPlease see help message by using 'perl 4b_minusInput.pl  -help' \n";
    exit 0;
}

## Get Arguments 
if ( exists $args{'-version' }   )     { say  "\n$version\n";    exit 0; }
if ( exists $args{'-help'    }   )     { say  "\n$HELP\n";       exit 0; }
if ( exists $args{'-in'      }   )     { $input_g  = $args{'-in'      }; }else{say   "\n -in   is required.\n";   say  "\n$HELP\n";    exit 0; }
if ( exists $args{'-out'     }   )     { $output_g = $args{'-out'     }; }else{say   "\n -out  is required.\n";   say  "\n$HELP\n";    exit 0; }

## Conditions
$input_g  =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";
$output_g =~ m/^\S+$/    ||  die   "\n\n$HELP\n\n";

## Print Command Arguments to Standard Output
say  "\n
        ################ Arguments ###############################
                Input   Path:  $input_g
                Output  Path:  $output_g
        ###############################################################
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

opendir(my $DH_input_g, "$input_g/IP")  ||  die;
my @inputFiles_g = readdir($DH_input_g);
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";




for ( my $i=0; $i<=$#inputFiles_g; $i++ ) {
        my $temp = $inputFiles_g[$i];
        next unless $temp =~ m/\.IP$/;
        next unless $temp =~ m/^sample/;
        my $IPfile = `ls   $input_g/IP/$temp/*.gene.rpkm.bed`;        
        $IPfile =~ s/\s//g or die;
        print "##$IPfile##\n";
        open(IP1, "<", $IPfile) or die;
        my @lines1 = <IP1>; 
        open(OUT1, ">", "$output_g/$temp.bed") or die;
        print OUT1 $lines1[0]; 
        
        $temp =~ s/\.IP/.IN/ or die;
        my $INPUTfile = `ls   $input_g/Input/$temp/*.gene.rpkm.bed`;        
        $INPUTfile =~ s/\s//g or die;
        print "##$INPUTfile##\n\n";
        open(IN1, "<", $INPUTfile) or die;
        my @lines2 = <IN1>; 
         
        ($#lines1 == $#lines2) or die;
        for (my $j=1; $j<=$#lines1; $j++) {
            $lines1[$j] =~ m/^(\d+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([-+])\t(\S+)\n$/   or die  "$lines1[$j]\n\n";
            my $chr1    = $1;
            my $start1  = $2;
            my $end1    = $3;
            my $name1   = $4;
            my $log1    = $5;
            my $strand1 = $6;
            my $fpkm1   = $7;
        
            $lines2[$j] =~ m/^(\d+)\t(\d+)\t(\d+)\t(\S+)\t(\S+)\t([-+])\t(\S+)\n$/   or die  "$lines2[$j]\n\n";
            my $chr2    = $1;
            my $start2  = $2;
            my $end2    = $3;
            my $name2   = $4;
            my $log2    = $5;
            my $strand2 = $6;
            my $fpkm2   = $7;
            
            $chr1 = "chr".$chr1;
            $chr2 = "chr".$chr2;
            my $score = $fpkm1 - $fpkm2; 
            if($score < 0) {$score = 0; }
            
            ( ($chr1 eq $chr2) and ($start1 == $start2) and ($end1 == $end2) and ($name1 eq $name2)  and ($strand1 eq $strand2) ) or die;
            print   OUT1   "$chr1\t$start1\t$end1\t$name1\t$log1\t$strand1\t$score\n";
        }
        
}



print("\n\n\n\n\n#########################################\n");
###################################################################################################################################################################################################






###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END




