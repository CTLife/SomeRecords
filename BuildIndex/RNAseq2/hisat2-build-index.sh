mkdir		hg38						
hisat2-build	-p 8    Shortcuts/hg38/hg38.fasta 	hg38/hg38	        >>hg38.buildIndex.runLog	2>&1

mkdir		ce11						
hisat2-build	-p 8    Shortcuts/ce11/ce11.fasta 	ce11/ce11	        >>ce11.buildIndex.runLog	2>&1

mkdir		danRer10						
hisat2-build	-p 8    Shortcuts/danRer10/danRer10.fasta 	danRer10/danRer10	        >>danRer10.buildIndex.runLog	2>&1

mkdir		dm6						
hisat2-build	-p 8    Shortcuts/dm6/dm6.fasta 	dm6/dm6	        >>dm6.buildIndex.runLog	2>&1


mkdir		mm10						
hisat2-build	-p 8    Shortcuts/mm10/mm10.fasta 	mm10/mm10	        >>mm10.buildIndex.runLog	2>&1


mkdir		sacCer3						
hisat2-build	-p 8    Shortcuts/sacCer3/sacCer3.fasta 	sacCer3/sacCer3	        >>sacCer3.buildIndex.runLog	2>&1





mkdir		OtherGenomes/Adapters						
hisat2-build	-p 8    Shortcuts/OtherGenomes/Adapters.fasta 	OtherGenomes/Adapters/Adapters	        >>Adapters.buildIndex.runLog	2>&1


mkdir		OtherGenomes/Ecoli						
hisat2-build	-p 8    Shortcuts/OtherGenomes/Ecoli.fasta 	OtherGenomes/Ecoli/Ecoli	        >>Ecoli.buildIndex.runLog	2>&1


mkdir		OtherGenomes/phiX174						
hisat2-build	-p 8    Shortcuts/OtherGenomes/phiX174.fasta 	OtherGenomes/phiX174/phiX174	        >>phiX174.buildIndex.runLog	2>&1


mkdir		OtherGenomes/UniVectors						
hisat2-build	-p 8    Shortcuts/OtherGenomes/UniVectors.fasta 	OtherGenomes/UniVectors/UniVectors	        >>UniVectors.buildIndex.runLog	2>&1




mkdir		hg38.ensembl						
hisat2-build	-p 8    Shortcuts/hg38/hg38.ensembl.genome.fasta 	hg38.ensembl/hg38.ensembl	        >>hg38.ensembl.buildIndex.runLog	2>&1

mkdir		mm10.ensembl						
hisat2-build	-p 8    Shortcuts/mm10/mm10.ensembl.genome.fasta 	mm10.ensembl/mm10.ensembl	        >>mm10.ensembl.buildIndex.runLog	2>&1



