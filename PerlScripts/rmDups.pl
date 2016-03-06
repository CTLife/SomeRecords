## remove duplicates based on the first column.



#!/usr/bin/env perl
use strict;
use warnings;




########## Keys and Values ##########
my %args = @ARGV;


########## Set Defaults ##########
my $input_g      = 'eBayes-down.xls';


########## Get Arguments ##########
if ( exists $args{'-i'}   )    { $input_g  = $args{'-i'};   }


########### Conditions #############
$input_g      =~ m/^\S+$/ or die;


######### Example ###########
# perl   rmDups.pl  -i  eBayes-down.xls   





my $output1_g = "z-noDups-$input_g";
my $output2_g = "z-Dups-$input_g";
open(INPUT_FH,      "<",   $input_g  )      or die "$!"; 
open(OUTPUT1_FH,    ">",   $output1_g)      or die "$!"; 
open(OUTPUT2_FH,    ">",   $output2_g)      or die "$!"; 


my @lines1 = <INPUT_FH>; 


for (my $i=0; $i<=$#lines1; $i++) {
    my $temp1 = $lines1[$i];
    $temp1 =~ m/^(\S+)\s+(\S+)\s+(\S+)\s*/ or die;
    my $id1 = $1;  ## based on this column
    my $bool = 0;
    for (my $j=$i+1; $j<=$#lines1; $j++) {
        my $temp2 = $lines1[$j];
        $temp2 =~ m/^(\S+)\s+(\S+)\s+(\S+)\s*/ or die;
        my $id2 = $1;  ## based on this column
        if($id1 eq $id2) {$bool = 1; }
    }
    if($bool == 0) {print  OUTPUT1_FH  $temp1;}
    if($bool == 1) {print  OUTPUT2_FH  $temp1;}
}




