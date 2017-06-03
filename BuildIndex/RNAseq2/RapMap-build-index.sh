mkdir		ce11			
rapmap  quasiindex   --transcripts   Shortcuts/ce11/ce11.RefSeq.fasta     --klen 31      --index ce11/ce11.RefSeq      --numThreads 4	     >>ce11.buildIndex.runLog 2>&1      


mkdir		danRer10			
rapmap  quasiindex   --transcripts   Shortcuts/danRer10/danRer10.RefSeq.fasta     --klen 31      --index danRer10/danRer10.RefSeq      --numThreads 4	     >>danRer10.buildIndex.runLog 2>&1      


mkdir		dm6			
rapmap  quasiindex   --transcripts   Shortcuts/dm6/dm6.RefSeq.fasta     --klen 31      --index dm6/dm6.RefSeq      --numThreads 4	     >>dm6.buildIndex.runLog 2>&1      


mkdir		hg38			
rapmap  quasiindex   --transcripts   Shortcuts/hg38/hg38.RefSeq.fasta     --klen 31      --index hg38/hg38.RefSeq      --numThreads 4	     >>hg38.buildIndex.runLog 2>&1      

mkdir		mm10			
rapmap  quasiindex   --transcripts   Shortcuts/mm10/mm10.RefSeq.fasta     --klen 31      --index mm10/mm10.RefSeq      --numThreads 4	     >>mm10.buildIndex.runLog 2>&1      


mkdir		sacCer3			
rapmap  quasiindex   --transcripts   Shortcuts/sacCer3/sacCer3.RefSeq.fasta     --klen 31      --index sacCer3/sacCer3.RefSeq      --numThreads 4	     >>sacCer3.buildIndex.runLog 2>&1      





mkdir		hg38.ensembl.cDNA			
rapmap  quasiindex   --transcripts   Shortcuts/hg38/hg38.ensembl.cDNA.fasta     --klen 31      --index hg38.ensembl.cDNA/hg38.ensembl.cDNA      --numThreads 4	     >>hg38.ensembl.cDNA.buildIndex.runLog 2>&1      

mkdir		mm10.ensembl.cDNA			
rapmap  quasiindex   --transcripts   Shortcuts/mm10/mm10.ensembl.cDNA.fasta     --klen 31      --index mm10.ensembl.cDNA/mm10.ensembl.cDNA      --numThreads 4	     >>mm10.ensembl.cDNA.buildIndex.runLog 2>&1      



