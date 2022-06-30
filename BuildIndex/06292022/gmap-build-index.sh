 


mkdir -p  "wuhCor1"                      
gmap_build      --dir=wuhCor1	 --genomedb=wuhCor1   ../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa        >>  wuhCor1.runLog  2>&1
  
mkdir -p  "danRer11"                      
gmap_build      --dir=danRer11	 --genomedb=danRer11   ../../../14_Genomes/UCSC/Zebrafish/danRer11.fa        >>  danRer11.runLog  2>&1
  
mkdir -p  "micMur2"                      
gmap_build      --dir=micMur2	 --genomedb=micMur2   ../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa        >>  micMur2.runLog  2>&1
  
mkdir -p  "chlSab2"                      
gmap_build      --dir=chlSab2	 --genomedb=chlSab2   ../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa        >>  chlSab2.runLog  2>&1
  
mkdir -p  "lambda_phage"                      
gmap_build      --dir=lambda_phage	 --genomedb=lambda_phage   ../../../14_Genomes/NCBI/lambda_phage/lambda_phage.fa        >>  lambda_phage.runLog  2>&1  
  
mkdir -p  "hg38"                                   
gmap_build      --dir=hg38	 --genomedb=hg38               ../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    >>  hg38.runLog  2>&1    

mkdir -p  "mm39"                                 
gmap_build      --dir=mm39	 --genomedb=mm39               ../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa    >>  mm39.runLog  2>&1  
         




