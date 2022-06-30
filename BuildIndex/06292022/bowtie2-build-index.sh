 
#################
mkdir -p  "wuhCor1"                      
bowtie2-build    ../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa      wuhCor1/wuhCor1      >>  wuhCor1.runLog  2>&1
  
mkdir -p  "danRer11"                      
bowtie2-build      ../../../14_Genomes/UCSC/Zebrafish/danRer11.fa    danRer11/danRer11      >>  danRer11.runLog  2>&1
  
mkdir -p  "micMur2"                      
bowtie2-build      ../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa    micMur2/micMur2      >>  micMur2.runLog  2>&1
  
mkdir -p  "chlSab2"                      
bowtie2-build     ../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa      chlSab2/chlSab2     >>  chlSab2.runLog  2>&1
  
mkdir -p  "lambda_phage"                      
bowtie2-build       ../../../14_Genomes/NCBI/lambda_phage/lambda_phage.fa   lambda_phage/lambda_phage      >>  lambda_phage.runLog  2>&1  
  
mkdir -p  "hg38"                                   
bowtie2-build       ../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    hg38/hg38    >>  hg38.runLog  2>&1    

mkdir -p  "mm39"                                 
bowtie2-build          ../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa  mm39/mm39    >>  mm39.runLog  2>&1  
         


