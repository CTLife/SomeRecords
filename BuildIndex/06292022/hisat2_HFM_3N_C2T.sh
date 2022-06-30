 
 
#################
mkdir -p  "wuhCor1"                      
hisat-3n-build	-p 8    --base-change C,T   ../../../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa      wuhCor1/wuhCor1      >>  wuhCor1.runLog  2>&1
  
mkdir -p  "danRer11"                      
hisat-3n-build	-p 8      --base-change C,T    ../../../../../14_Genomes/UCSC/Zebrafish/danRer11.fa    danRer11/danRer11      >>  danRer11.runLog  2>&1
  
mkdir -p  "micMur2"                      
hisat-3n-build	-p 8      --base-change C,T    ../../../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa    micMur2/micMur2      >>  micMur2.runLog  2>&1
  
mkdir -p  "chlSab2"                      
hisat-3n-build	-p 8      --base-change C,T   ../../../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa      chlSab2/chlSab2     >>  chlSab2.runLog  2>&1
  
mkdir -p  "lambda_phage"                      
hisat-3n-build	-p 8      --base-change C,T     ../../../../../14_Genomes/NCBI/lambda_phage/lambda_phage.fa   lambda_phage/lambda_phage      >>  lambda_phage.runLog  2>&1  
  
mkdir -p  "hg38"                                   
hisat-3n-build	-p 8       --base-change C,T    ../../../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    hg38/hg38    >>  hg38.runLog  2>&1    

mkdir -p  "mm39"                                 
hisat-3n-build	-p 8       --base-change C,T     ../../../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa  mm39/mm39    >>  mm39.runLog  2>&1  
         

 

 
 
#################
mkdir -p  "wuhCor1_repeat"                      
hisat-3n-build	-p 8      --base-change C,T      --repeat-index     ../../../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa      wuhCor1_repeat/wuhCor1_repeat      >>  wuhCor1_repeat.runLog  2>&1
  
mkdir -p  "danRer11_repeat"                      
hisat-3n-build	-p 8      --base-change C,T       --repeat-index     ../../../../../14_Genomes/UCSC/Zebrafish/danRer11.fa    danRer11_repeat/danRer11_repeat      >>  danRer11_repeat.runLog  2>&1
  
mkdir -p  "micMur2_repeat"                      
hisat-3n-build	-p 8      --base-change C,T       --repeat-index     ../../../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa    micMur2_repeat/micMur2_repeat      >>  micMur2_repeat.runLog  2>&1
  
mkdir -p  "chlSab2_repeat"                      
hisat-3n-build	-p 8      --base-change C,T      --repeat-index     ../../../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa      chlSab2_repeat/chlSab2_repeat     >>  chlSab2_repeat.runLog  2>&1
  
mkdir -p  "lambda_phage_repeat"                      
hisat-3n-build	-p 8      --base-change C,T      --repeat-index     ../../../../../14_Genomes/NCBI/lambda_phage/lambda_phage.fa   lambda_phage_repeat/lambda_phage_repeat      >>  lambda_phage_repeat.runLog  2>&1  
  
mkdir -p  "hg38_repeat"                                   
hisat-3n-build	-p 8       --base-change C,T     --repeat-index     ../../../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    hg38_repeat/hg38_repeat    >>  hg38_repeat.runLog  2>&1    

mkdir -p  "mm39_repeat"                                 
hisat-3n-build	-p 8       --base-change C,T     --repeat-index     ../../../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa  mm39_repeat/mm39_repeat    >>  mm39_repeat.runLog  2>&1  
         

 
