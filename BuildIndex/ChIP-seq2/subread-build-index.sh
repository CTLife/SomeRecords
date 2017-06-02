mkdir  mm10
subread-buildindex  -o 	mm10/mm10   Shortcuts/mm10/mm10.fasta   >>mm10.buildIndex.runLog 2>&1                   

mkdir  hg38
subread-buildindex  -o 	hg38/hg38   Shortcuts/hg38/hg38.fasta   >>hg38.buildIndex.runLog 2>&1                   

mkdir  sacCer3
subread-buildindex  -o 	sacCer3/sacCer3   Shortcuts/sacCer3/sacCer3.fasta   >>sacCer3.buildIndex.runLog 2>&1                   


mkdir  ce11
subread-buildindex  -o 	ce11/ce11   Shortcuts/ce11/ce11.fasta   >>ce11.buildIndex.runLog 2>&1                   

mkdir  danRer10
subread-buildindex  -o 	danRer10/danRer10   Shortcuts/danRer10/danRer10.fasta   >>danRer10.buildIndex.runLog 2>&1                   

mkdir  dm6
subread-buildindex  -o 	dm6/dm6   Shortcuts/dm6/dm6.fasta   >>dm6.buildIndex.runLog 2>&1                   




mkdir -p  OtherGenomes/Adapters
subread-buildindex  -o 	OtherGenomes/Adapters/Adapters   Shortcuts/OtherGenomes/Adapters.fasta   >>Adapters.buildIndex.runLog 2>&1                   

mkdir -p  OtherGenomes/Ecoli
subread-buildindex  -o 	OtherGenomes/Ecoli/Ecoli   Shortcuts/OtherGenomes/Ecoli.fasta   >>Ecoli.buildIndex.runLog 2>&1                   

mkdir -p  OtherGenomes/phiX174
subread-buildindex  -o 	OtherGenomes/phiX174/phiX174   Shortcuts/OtherGenomes/phiX174.fasta   >>phiX174.buildIndex.runLog 2>&1                   

mkdir -p  OtherGenomes/UniVectors
subread-buildindex  -o 	OtherGenomes/UniVectors/UniVectors   Shortcuts/OtherGenomes/UniVectors.fasta   >>UniVectors.buildIndex.runLog 2>&1                   




mkdir  mm10.ensembl
subread-buildindex  -o 	mm10.ensembl/mm10.ensembl   Shortcuts/mm10/mm10.ensembl.genome.fasta   >>mm10.ensembl.buildIndex.runLog 2>&1                   

mkdir  hg38.ensembl
subread-buildindex  -o 	hg38.ensembl/hg38.ensembl   Shortcuts/hg38/hg38.ensembl.genome.fasta   >>hg38.ensembl.buildIndex.runLog 2>&1                   



