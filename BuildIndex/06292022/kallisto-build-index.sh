


mkdir -p  "danRer11"                      
kallisto  index      --index=danRer11/danRer11   ../../../14_Genomes/UCSC/Zebrafish/mrna.fa        >>  danRer11.runLog  2>&1
  
mkdir -p  "micMur2"                      
kallisto  index      --index=micMur2/micMur2   ../../../14_Genomes/UCSC/Mouse_lemur/mrna.fa        >>  micMur2.runLog  2>&1
  
mkdir -p  "chlSab2"                      
kallisto  index      --index=chlSab2/chlSab2   ../../../14_Genomes/UCSC/Green_monkey/mrna.fa        >>  chlSab2.runLog  2>&1
  
 
mkdir -p  "hg38"                                   
kallisto  index      --index=hg38/hg38               ../../../14_Genomes/Gencode/Human/gencode.v40.transcripts.fa    >>  hg38.runLog  2>&1    

mkdir -p  "mm10"                                 
kallisto  index      --index=mm10/mm10               ../../../14_Genomes/Gencode/Mouse/gencode.vM29.transcripts.fa    >>  mm10.runLog  2>&1  
         


