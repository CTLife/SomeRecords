#!/usr/bin/env  perl5
use  strict;
use  warnings;
use  v5.16;
## Perl5 version >= 5.16
## You can create a symbolic link for perl5 by using "sudo  ln  /usr/bin/perl   /usr/bin/perl5" in Ubuntu.
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

my $input_g = "1-raw";
my $output_g = "2-format";
my $output2_g = "3-for-Figure";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}
&myMakeDir($output_g);
&myMakeDir($output2_g);

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
        open(OUTPUT2_FH,  ">",   "$output2_g/$temp" )      or   die "$!"; 
        my @lines1 = <INPUT1_FH>; 

        my $Promoter = 0;
        my $Exon = 0;
        my $Intron = 0;
        my $Intergenic = 0;
        my $Others = 0;

        for (my $i=4; $i<=$#lines1; $i++) {
            my $temp1 = $lines1[$i];
            $temp1 =~ m/^(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s*/ or die;
            my $name1 = $2;
            my $name2 = $3;
            my $Frequency = $4;
            print  OUTPUT1_FH   "$name1.$name2\t$Frequency\n";
            if($name1 =~ m/Promoter/)   {$Promoter = $Promoter + $Frequency; }
            if($name2 =~ m/UTR/)        {$Others   = $Others + $Frequency; }
            if($name2 =~ m/Exon/)       {$Exon     = $Exon + $Frequency; }
            if($name2 =~ m/Intron/)     {$Intron   = $Intron + $Frequency; }
            if($name1 =~ m/Downstream/) {$Others   = $Others + $Frequency; }
            if($name1 =~ m/Distal/)     {$Intergenic = $Intergenic + $Frequency; }
        }
     
        print  OUTPUT2_FH   "Promoter\t$Promoter\n";
        print  OUTPUT2_FH   "Exon\t$Exon\n";
        print  OUTPUT2_FH   "Intron\t$Intron\n";
        print  OUTPUT2_FH   "Intergenic\t$Intergenic\n";
        print  OUTPUT2_FH   "Others\t$Others\n";

}


}
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "\tJob Done! Cheers! \n\n\n\n\n";





## END
