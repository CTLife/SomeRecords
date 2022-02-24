input1=6_pca
out1=8_cis_permutation
mkdir  $out1
file1=merged

 
QTLtools cis --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge/merged.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --permute 1000    --out $out1/$file1.nominals.txt  --log $out1/$file1.Log   >  $out1/$file1.runLog.txt 2>&1 
 


QTLtools cis --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge/merged.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --permute 10000    --out $out1/$file1.2.nominals.txt  --log $out1/$file1.2.Log   >  $out1/$file1.2.runLog.txt 2>&1 
 

 
