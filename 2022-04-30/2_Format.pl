#!/usr/bin/env  perl
use  strict;
use  warnings;
use  v5.12;
## Perl5 version >= 5.12
## Suffixes of all self-defined global variables must be "_g".
###################################################################################################################################################################################################



my $input1="1_raw_m6A-QTLs";
my $out="2_Format";


###################################################################################################################################################################################################
say   "\n\n\n\n\n\n##################################################################################################";
say   "Running ......";

sub myMakeDir  {
    my $path = $_[0];
    if ( !( -e $path) )  { system("mkdir  -p  $path"); }
    if ( !( -e $path) )  { mkdir $path  ||  die;       }
}
&myMakeDir($out);

opendir(my $FH_input1_g, $input1)  ||  die;
my @files = readdir($FH_input1_g);
###################################################################################################################################################################################################



 
        
###################################################################################################################################################################################################

for(my $i1=0; $i1<=$#files; $i1=$i1+1) {
        my $temp1 = $files[$i1];
        next unless $temp1 =~ m/\.txt$/;
        my $file1 = "$input1/$temp1";
        my $file2 = "$out/$temp1";

        open(my $FH1, "<", $file1)  ||  die;
        open(my $FH2, ">", $file2)  ||  die;
        print  $FH2   "phe_id\tnew_var_id\tphe_chr\tphe_from\tphe_to\tphe_strd\tn_var_in_cis\tdist_phe_var\tvar_id\tvar_chr\tvar_from\tvar_to\tnom_pval\tr_squared\tslope\tslope_se\tbest_hit\n";
         
        while(my $line=<$FH1>) {
           my $temp2=$line; 
           $temp2 =~ m/^(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+([-+])\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+([01])\n$/  or  die "#$temp2#\n"; 
           my $phe_id        = $1;
           my $phe_chr       = $2;
           my $phe_from      = $3;
           my $phe_to        = $4;
           my $phe_strd      = $5;
           my $n_var_in_cis  = $6;
           my $dist_phe_var  = $7;
           my $var_id        = $8;	
           my $var_chr       = $9;
           my $var_from      = $10;
           my $var_to        = $11;
           my $nom_pval      = $12;
           my $r_squared     = $13;
           my $slope         = $14;
           my $slope_se      = $15;
           my $best_hit      = $16;
           my $new_var_id =  "SNP"."_". "$var_id"."_"."$var_chr"."_"."$var_from"."_"."$var_to" ;
           $new_var_id =~ s/_\._/_NA_/;
           print  $FH2   "$phe_id\t$new_var_id\t$phe_chr\t$phe_from\t$phe_to\t$phe_strd\t$n_var_in_cis\t$dist_phe_var\t$var_id\t$var_chr\t$var_from\t$var_to\t$nom_pval\t$r_squared\t$slope\t$slope_se\t$best_hit\n";
        }
}

###################################################################################################################################################################################################






say   "Done ......";






#####
