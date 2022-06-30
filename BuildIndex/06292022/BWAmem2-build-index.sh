 
#################
mkdir -p  "wuhCor1"                      
bwa-mem2  index      -p wuhCor1/wuhCor1   ../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa        >>  wuhCor1.runLog  2>&1
  
mkdir -p  "danRer11"                      
bwa-mem2  index      -p danRer11/danRer11   ../../../14_Genomes/UCSC/Zebrafish/danRer11.fa        >>  danRer11.runLog  2>&1
  
mkdir -p  "micMur2"                      
bwa-mem2  index      -p micMur2/micMur2   ../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa        >>  micMur2.runLog  2>&1
  
mkdir -p  "chlSab2"                      
bwa-mem2  index      -p chlSab2/chlSab2   ../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa        >>  chlSab2.runLog  2>&1
  
mkdir -p  "lambda_phage"                      
bwa-mem2  index      -p lambda_phage/lambda_phage   ../../../14_Genomes/NCBI/lambda_phage/lambda_phage.fa        >>  lambda_phage.runLog  2>&1  
  
mkdir -p  "hg38"                                   
bwa-mem2  index      -p hg38/hg38               ../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    >>  hg38.runLog  2>&1    

mkdir -p  "mm39"                                 
bwa-mem2  index      -p mm39/mm39               ../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa    >>  mm39.runLog  2>&1  
         


