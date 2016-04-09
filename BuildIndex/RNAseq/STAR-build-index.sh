mkdir		ce11			
mkdir		mm10			
												
STAR  --runThreadN 4   --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/ce11.fa		--genomeDir  ce11	--sjdbGTFfile  Shortcuts/ce11_RefSeq_GTF  	--sjdbOverhang 99	>>ce11.buildIndex.runLog	2>&1          
STAR  --runThreadN 4   --runMode genomeGenerate  	--genomeFastaFiles   Shortcuts/mm10.fa		--genomeDir  mm10	--sjdbGTFfile  Shortcuts/mm10_RefSeq_GTF  	--sjdbOverhang 99	>>mm10.buildIndex.runLog	2>&1            



