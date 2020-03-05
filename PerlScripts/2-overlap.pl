#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################





###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}


my $output_g = "2-overlap";
my $output2_g = "3-finalPeaks/byIDR";
my $input_g = "1-rawData";
&myMakeDir($output_g);
&myMakeDir($output2_g);

opendir(FH1, $input_g)  ||  die;
###################################################################################################################################################################################################





###################################################################################################################################################################################################
my @folders = readdir(FH1);

for(my $i=0; $i<=$#folders; $i++) {
    my $folder = $folders[$i]; 
    next unless $folder !~ m/^[.]/;
    next unless $folder !~ m/[~]$/;
    next unless $folder !~ m/\.sh$/;
    ##say "$input_g/$folder";
    opendir(FH2, "$input_g/$folder")  ||  die;
    my @files = readdir(FH2);
    @files = grep(/\.bed$/, @files);
    ##@files = grep -f, @files;
    print join("\t",@files),"\n";
    
    ($#files == 2) or die "##@files##\n\n";
    my @sortfiles = sort{-s "$input_g/$folder/$a" <=> -s "$input_g/$folder/$b"} @files; 
    ($#sortfiles == 2) or die;
    system("cp     $input_g/$folder/$sortfiles[0]      $output2_g/$folder.bed"  );
    system(" intervene  venn  --input  $input_g/$folder/$sortfiles[1]  $input_g/$folder/$sortfiles[2]  --names=rep1,rep2   --bedtools-options s   --output $output_g/$folder      --save-overlaps      --dpi 1200    --figtype svg ");                  
    print "########################################\n\n\n\n\n";

}


 







#####
