#!/usr/bin/env perl
use strict;
use warnings;





########## Keys and Values ##########
my %args = @ARGV;


########## Set Defaults ##########
my $file1_g      = 'GSE53386_Mouse_RNA.txt';
my $file2_g      = 'GSE79552.FPKM.mouse.txt';


########## Get Arguments ##########
if ( exists $args{'-file1'}   )    { $file1_g      = $args{'-file1'};   }
if ( exists $args{'-file2'}   )    { $file2_g      = $args{'-file2'};   }


########### Conditions #############
$file1_g      =~ m/^\S+$/ or die;
$file2_g      =~ m/^\S+$/ or die;


######### Example ###########
# perl   compare.pl   -file1  GSE53386_Mouse_RNA.txt     -file2   GSE79552.FPKM.mouse.txt  
####################################################################################################





####################################################################################################
open(INPUT1,    "<",   "$file1_g")     or die "$!"; 
open(INPUT2,    "<",   "$file2_g")     or die "$!"; 
my @lines1 = <INPUT1>; 
my @lines2 = <INPUT2>; 

my $common1    = "z--common1--$file1_g--$file2_g";
my $common2    = "z--common2--$file2_g--$file1_g";
my $notFound1  = "z--notFound1--$file1_g";
my $notFound2  = "z--notFound2--$file2_g";
open(common1_FH,     ">",      $common1    )     or  die "$!";   
open(common2_FH,     ">",      $common2    )     or  die "$!";   
open(notFound1_FH,   ">",      $notFound1  )     or  die "$!";   
open(notFound2_FH,   ">",      $notFound2  )     or  die "$!";   


for (my $i=0; $i<=$#lines1; $i++) {
    my $bool = 0; 
    my $temp1 = $lines1[$i];
    $temp1 =~ m/^(\S+)\s/ or die;
    my $ID1 = $1;
    $temp1 =~ s/\n// ;
    print  "#$ID1#\n";
    for (my $j=0; $j<=$#lines2; $j++) {
        my $temp2 = $lines2[$j];
        $temp2 =~ s/\n// ;
        $temp2 =~ m/^(\S+)\s+/ or die "\n$temp2\n\n";
        my $id2 = $1; 
        if($ID1 =~ m/^$id2$/i)  {print  common1_FH  "$temp1\t$temp2\n"; $bool=1; last; }
    }
    if($bool == 0) { print    notFound1_FH   "$temp1\tz-noFound\n";}
}




for (my $i=0; $i<=$#lines2; $i++) {
    my $bool = 0; 
    my $temp1 = $lines2[$i];
    $temp1 =~ m/^(\S+)\s/ or die;
    my $ID1 = $1;
    $temp1 =~ s/\n// ;
    print  "#$ID1#\n";
    for (my $j=0; $j<=$#lines1; $j++) {
        my $temp2 = $lines1[$j];
        $temp2 =~ s/\n// ;
        $temp2 =~ m/^(\S+)\s+/ or die "\n$temp2\n\n";
        my $id2 = $1; 
        if($ID1 =~ m/^$id2$/i)  {print  common2_FH  "$temp1\t$temp2\n"; $bool=1; last; }
    }
    if($bool == 0) { print    notFound2_FH   "$temp1\tz-noFound\n";}
}



####################################################################################################




