input1=6_pca
out1=10_trans_full
mkdir  $out1
file1=merged

 
QTLtools trans --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge/merged.sorted.bed.gz       \
            --nominal --threshold 1e-2     --out $out1/$file1.nominals.txt  --log $out1/$file1.Log   >  $out1/$file1.runLog.txt 2>&1 

 
QTLtools trans --vcf   ReadMe/20220210/genotype.chr1-22.vcf.gz   --bed  5_merge/merged.sorted.bed.gz       \
            --permute --threshold 1e-2     --out $out1/$file1.trans.perm123.txt  --log $out1/$file1.trans.perm123.Log  --seed 123  >  $out1/$file1.runLog.txt 2>&1 
 
  



 
