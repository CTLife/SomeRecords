#!/usr/bin/env perl
use strict;
use warnings;





########## Keys and Values ##########
my %args = @ARGV;



########## Set Defaults ##########
my $file1_g      = '4-EEDKO-Down.xls';
my $file2_g      = '5-Down-GuFei.bed';


########## Get Arguments ##########
if ( exists $args{'-file1'}   )    { $file1_g      = $args{'-file1'};   }
if ( exists $args{'-file2'}   )    { $file2_g      = $args{'-file2'};   }


########### Conditions #############
$file1_g      =~ m/^\S+$/ or die;
$file2_g      =~ m/^\S+$/ or die;


######### Example ###########
# perl   findName1.pl   -file1  GenesOfTFs      -file2   allH3K27acPeaks.100kb_symbol  






open(INPUT1,    "<",   "$file1_g")     or die "$!"; 
open(INPUT2,    "<",   "$file2_g")     or die "$!"; 

my $common1File = "$file1_g-$file2_g";
my $common2File = "noFind-$file1_g-$file2_g";
open(Common1,   ">",   "$common1File")     or die "$!";   
open(Common2,   ">",   "$common2File")     or die "$!";   



my @lines1 = <INPUT1>; 
my @lines2 = <INPUT2>; 


for (my $i=0; $i<=$#lines1; $i++) {
    my $bool = 0; 
    my $temp1 = $lines1[$i];
    $temp1 =~ m/^([^\t]+)\t/ or die;
    my $ID1 = $1;
    $temp1 =~ s/\n// ;
    print  "#$ID1#\n";
    $ID1  = "\t".$ID1.",,";
    for (my $j=0; $j<=$#lines2; $j++) {
        my $temp2 = $lines2[$j];
        $temp2 =~ s/\n// ;
        $temp2 =~ m/^(\S+)\s+/ or die "\n$temp2\n\n";
        my $id2 = $1; 
        if($ID1 =~ m/\s$id2,/i)  {print  Common1  "$temp1\t$temp2\n"; $bool=1;  }
    }
    if($bool == 0) { print    Common2   "$temp1\tz-noFound\n";}

}









