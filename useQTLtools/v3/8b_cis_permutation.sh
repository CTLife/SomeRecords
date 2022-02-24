input1=6b_pca
out1=8b_cis_permutation
mkdir  $out1
file1=merged

 
QTLtools cis --vcf   otherFiles/genotype.chr1-22.vcf.gz   --bed  5b_merge/$file1.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --permute 1000    --out $out1/$file1.permutation.txt  --log $out1/$file1.Log   >  $out1/$file1.runLog.txt 2>&1 
 


QTLtools cis --vcf   otherFiles/genotype.chr1-22.vcf.gz   --bed  5b_merge/$file1.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --permute 10000    --out $out1/$file1.2.permutation.txt  --log $out1/$file1.2.Log   >  $out1/$file1.2.runLog.txt 2>&1 
 

 
