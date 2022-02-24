input1=6_pca
out1=7_cis_nominal/minus
mkdir  $out1
file1=merged

 
QTLtools cis --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge_peaks/minus.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --nominal 1    --out $out1/$file1.nominals.txt  --log $out1/$file1.Log   >  $out1/$file1.runLog.txt 2>&1 


QTLtools cis --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge_peaks/minus.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --nominal 0.1    --out $out1/$file1.2.nominals.txt  --log $out1/$file1.2.Log   >  $out1/$file1.2.runLog.txt 2>&1 


QTLtools cis --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge_peaks/minus.sorted.bed.gz   --cov   $input1/$file1.sorted.bed.pca   \
            --nominal 0.01    --out $out1/$file1.3.nominals.txt  --log $out1/$file1.3.Log   >  $out1/$file1.3.runLog.txt 2>&1 





 
