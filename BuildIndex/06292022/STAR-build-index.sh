


mkdir		SARS-CoV-2														
STAR  --runThreadN 10    --runMode genomeGenerate  	--genomeFastaFiles   ../../../14_Genomes/UCSC/SARS-CoV-2/wuhCor1.fa	 --genomeDir  SARS-CoV-2	 >> SARS-CoV-2.buildIndex.runLog	2>&1            


mkdir		mm39														
STAR  --runThreadN 10    --runMode genomeGenerate  	--genomeFastaFiles  ../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa	 --genomeDir  mm39	--sjdbGTFfile  ../../../14_Genomes/Gencode/Mouse/gencode.vM29.basic.annotation.gtf 	>> mm39.buildIndex.runLog	2>&1            


mkdir		hg38														
STAR  --runThreadN 10    --runMode genomeGenerate  	--genomeFastaFiles   ../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa    --genomeDir  hg38	 --sjdbGTFfile  ../../../14_Genomes/Gencode/Human/gencode.v40.basic.annotation.gtf	>> hg38.buildIndex.runLog	2>&1                                                                              





mkdir		danRer11														
STAR  --runThreadN 10    --runMode genomeGenerate  	--genomeFastaFiles   ../../../14_Genomes/UCSC/Zebrafish/danRer11.fa	 --genomeDir  danRer11	  --sjdbGTFfile ../../../14_Genomes/UCSC/Zebrafish/genes/danRer11.ncbiRefSeq.gtf    >> danRer11.buildIndex.runLog	2>&1            




mkdir		micMur2														
STAR  --runThreadN 10    --runMode genomeGenerate  	--genomeFastaFiles   ../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa	 --genomeDir  micMur2	  --sjdbGTFfile ../../../14_Genomes/UCSC/Mouse_lemur/genes/micMur2.ncbiRefSeq.gtf    >> micMur2.buildIndex.runLog	2>&1            




mkdir		chlSab2														
STAR  --runThreadN 10    --runMode genomeGenerate  	--genomeFastaFiles   ../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa	 --genomeDir  chlSab2	  --sjdbGTFfile ../../../14_Genomes/UCSC/Green_monkey/genes/chlSab2.ncbiRefSeq.gtf    >> chlSab2.buildIndex.runLog	2>&1            




