mkdir		ce11			
salmon index   --transcripts   Shortcuts/ce11/ce11.RefSeq.fasta     --kmerLen 31      --index ce11/ce11.RefSeq      --threads  8    	  --type quasi    	     >>ce11.buildIndex.runLog 2>&1      


mkdir		danRer10			
salmon index   --transcripts   Shortcuts/danRer10/danRer10.RefSeq.fasta     --kmerLen 31      --index danRer10/danRer10.RefSeq      --threads  8    	  --type quasi    	     >>danRer10.buildIndex.runLog 2>&1      


mkdir		dm6			
salmon index   --transcripts   Shortcuts/dm6/dm6.RefSeq.fasta     --kmerLen 31      --index dm6/dm6.RefSeq      --threads  8    	  --type quasi    	     >>dm6.buildIndex.runLog 2>&1      


mkdir		hg38			
salmon index   --transcripts   Shortcuts/hg38/hg38.RefSeq.fasta     --kmerLen 31      --index hg38/hg38.RefSeq      --threads  8    	  --type quasi    	     >>hg38.buildIndex.runLog 2>&1      

mkdir		mm10			
salmon index   --transcripts   Shortcuts/mm10/mm10.RefSeq.fasta     --kmerLen 31      --index mm10/mm10.RefSeq      --threads  8    	  --type quasi    	     >>mm10.buildIndex.runLog 2>&1      


mkdir		sacCer3			
salmon index   --transcripts   Shortcuts/sacCer3/sacCer3.RefSeq.fasta     --kmerLen 31      --index sacCer3/sacCer3.RefSeq      --threads  8    	  --type quasi    	     >>sacCer3.buildIndex.runLog 2>&1      





mkdir		hg38.ensembl.cDNA			
salmon index   --transcripts   Shortcuts/hg38/hg38.ensembl.cDNA.fasta     --kmerLen 31      --index hg38.ensembl.cDNA/hg38.ensembl.cDNA      --threads  8    	  --type quasi    	     >>hg38.ensembl.cDNA.buildIndex.runLog 2>&1      

mkdir		mm10.ensembl.cDNA			
salmon index   --transcripts   Shortcuts/mm10/mm10.ensembl.cDNA.fasta     --kmerLen 31      --index mm10.ensembl.cDNA/mm10.ensembl.cDNA      --threads  8    	  --type quasi    	     >>mm10.ensembl.cDNA.buildIndex.runLog 2>&1      



