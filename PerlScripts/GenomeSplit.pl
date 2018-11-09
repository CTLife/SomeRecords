#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;   ## for using fileparse



########## Keys and Values ##########
my %args = @ARGV;


########## Set Defaults ##########
my $input_g = 'hg38.chrom.sizes';
my $bin_g   = 1000;

########## Get Arguments ##########
if ( exists $args{'-i'}   )    { $input_g  = $args{'-i'};   }
if ( exists $args{'-s'}   )    { $bin_g    = $args{'-s'};   }


########### Conditions #############
$input_g    =~ m/^\S+$/ or die;
$bin_g      =~ m/^\d+$/ or die;


######### Example ###########
# perl   GenomeSplit.pl     -i  hg38.chrom.sizes     -s  1000  




my $output1_g = $input_g.".$bin_g"."bpTiles.bed";
open(OUTPUT_FH, ">",   $output1_g  )        or die "$!"; 
open(INPUT_FH,  "<",   $input_g    )        or die "$!"; 
my @lines1 = <INPUT_FH>; 


my  $num = 0;
for (my $i=0; $i<=$#lines1; $i++) {
    my $temp1 = $lines1[$i];
    $temp1 =~ m/^(\S+)\s+(\d+)\s*/ or die;
    my $chr    = $1;   
    my $total  = $2;
    for (my $j=1; $j<=$total; $j=$j+$bin_g) {
        my $start = $j;
        my $end   = $j+$bin_g;
        if($end > $total) {$end = $total;} 
        $num = $num + 1;
        my $name = "Tiles".".$bin_g"."bp_$num"; 
        print  OUTPUT_FH  "$chr\t$start\t$end\t$name\n";
    }

}







