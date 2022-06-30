 
####################
mkdir -p  "wuhCor1"                      
subread-buildindex    -F   -o wuhCor1/wuhCor1   ../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa        >>  wuhCor1.runLog  2>&1
  
mkdir -p  "danRer11"                      
subread-buildindex    -F   -o danRer11/danRer11   ../../../14_Genomes/UCSC/Zebrafish/danRer11.fa        >>  danRer11.runLog  2>&1
  
mkdir -p  "micMur2"                      
subread-buildindex   -F    -o micMur2/micMur2   ../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa        >>  micMur2.runLog  2>&1
  
mkdir -p  "chlSab2"                      
subread-buildindex   -F    -o chlSab2/chlSab2   ../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa        >>  chlSab2.runLog  2>&1
  
mkdir -p  "lambda_phage"                      
subread-buildindex   -F    -o lambda_phage/lambda_phage   ../../../14_Genomes/NCBI/lambda_phage/lambda_phage.fa        >>  lambda_phage.runLog  2>&1  
  
mkdir -p  "hg38"                                   
subread-buildindex   -F    -o hg38/hg38               ../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    >>  hg38.runLog  2>&1    

mkdir -p  "mm39"                                 
subread-buildindex    -F   -o mm39/mm39               ../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa    >>  mm39.runLog  2>&1  
         


