
mkdir		mm10														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/mm10/mm10.fasta		--genomeDir  mm10	   --sjdbGTFfile  Shortcuts/mm10/mm10.RefSeq.GTF    	--sjdbOverhang 99	>>mm10.buildIndex.runLog	2>&1            



mkdir		ce11														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/ce11/ce11.fasta		--genomeDir  ce11	   --sjdbGTFfile  Shortcuts/ce11/ce11.RefSeq.GTF    	--sjdbOverhang 99	>>ce11.buildIndex.runLog	2>&1            



mkdir		danRer10														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/danRer10/danRer10.fasta		--genomeDir  danRer10	   --sjdbGTFfile  Shortcuts/danRer10/danRer10.RefSeq.GTF    	--sjdbOverhang 99	>>danRer10.buildIndex.runLog	2>&1            



mkdir		dm6														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/dm6/dm6.fasta		--genomeDir  dm6	   --sjdbGTFfile  Shortcuts/dm6/dm6.RefSeq.GTF    	--sjdbOverhang 99	>>dm6.buildIndex.runLog	2>&1            



mkdir		hg38														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/hg38/hg38.fasta		--genomeDir  hg38	   --sjdbGTFfile  Shortcuts/hg38/hg38.RefSeq.GTF    	--sjdbOverhang 99	>>hg38.buildIndex.runLog	2>&1            



mkdir		sacCer3														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/sacCer3/sacCer3.fasta		--genomeDir  sacCer3	   --sjdbGTFfile  Shortcuts/sacCer3/sacCer3.RefSeq.GTF    	--sjdbOverhang 99	>>sacCer3.buildIndex.runLog	2>&1            










mkdir		mm10.ensembl														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/mm10/mm10.ensembl.genome.fasta		--genomeDir  mm10.ensembl	   --sjdbGTFfile  Shortcuts/mm10/mm10.ensembl.cDNA.GTF    	--sjdbOverhang 99	>>mm10.ensembl.buildIndex.runLog	2>&1            

mkdir		hg38.ensembl														
STAR  --runThreadN 20    --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/hg38/hg38.ensembl.genome.fasta		--genomeDir  hg38.ensembl	   --sjdbGTFfile  Shortcuts/hg38/hg38.ensembl.cDNA.GTF    	--sjdbOverhang 99	>>hg38.ensembl.buildIndex.runLog	2>&1            





