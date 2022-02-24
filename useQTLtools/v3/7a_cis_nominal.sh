input1=6a_pca
out1=7a_cis_nominal
mkdir  $out1
file1=merged

 
QTLtools cis --vcf   otherFiles/genotype.chr1-22.vcf.gz   --bed  5a_merge/$file1.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --nominal 0.01    --out $out1/$file1.nominals.txt  --log $out1/$file1.Log   >  $out1/$file1.runLog.txt 2>&1 







 
