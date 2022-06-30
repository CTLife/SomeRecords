
# Build the standard HISAT-3N index  
mkdir  hg38_C2T
hisat-3n-build  -p 8   --base-change C,T                      ../../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa   hg38_C2T/hg38 

# Build the repeat HISAT-3N index  
mkdir  hg38_C2T_repeat
hisat-3n-build  -p 8   --base-change C,T   --repeat-index     ../../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa   hg38_C2T_repeat/hg38  


