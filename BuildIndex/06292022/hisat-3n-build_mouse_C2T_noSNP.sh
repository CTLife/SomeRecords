# Build the standard HISAT-3N index 
mkdir  mm39_C2T
hisat-3n-build  -p 8   --base-change C,T       ../../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa   mm39_C2T/mm39 

# Build the repeat HISAT-3N index  
mkdir  mm39_C2T_repeat
hisat-3n-build  -p 8   --base-change C,T   --repeat-index     ../../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa   mm39_C2T_repeat/mm39  



