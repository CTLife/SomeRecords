mkdir  -p  sacCer3                					
lastdb -uBISF        sacCer3/sacCer3_forward        Shortcuts/sacCer3/sacCer3.fasta	      >>sacCer3.buildIndex.runLog	2>&1
lastdb -uBISR        sacCer3/sacCer3_reverse        Shortcuts/sacCer3/sacCer3.fasta	      >>sacCer3.buildIndex.runLog	2>&1

mkdir -p mm10
lastdb -uBISF        mm10/mm10_forward	        Shortcuts/mm10/mm10.fasta	                      >>mm10.buildIndex.runLog	2>&1
lastdb -uBISR        mm10/mm10_reverse	        Shortcuts/mm10/mm10.fasta	                      >>mm10.buildIndex.runLog	2>&1

mkdir -p hg38
lastdb -uBISF        hg38/hg38_forward 	        Shortcuts/hg38/hg38.fasta	                      >>hg38.buildIndex.runLog	2>&1
lastdb -uBISR        hg38/hg38_reverse 	        Shortcuts/hg38/hg38.fasta	                      >>hg38.buildIndex.runLog	2>&1

					
mkdir  -p  ce11                					
lastdb -uBISF        ce11/ce11_forward        Shortcuts/ce11/ce11.fasta	      >>ce11.buildIndex.runLog	2>&1
lastdb -uBISR        ce11/ce11_reverse        Shortcuts/ce11/ce11.fasta	      >>ce11.buildIndex.runLog	2>&1

mkdir -p danRer10
lastdb -uBISF        danRer10/danRer10_forward	        Shortcuts/danRer10/danRer10.fasta	                      >>danRer10.buildIndex.runLog	2>&1
lastdb -uBISR        danRer10/danRer10_reverse	        Shortcuts/danRer10/danRer10.fasta	                      >>danRer10.buildIndex.runLog	2>&1

mkdir -p dm6
lastdb -uBISF        dm6/dm6_forward 	        Shortcuts/dm6/dm6.fasta	                      >>dm6.buildIndex.runLog	2>&1
lastdb -uBISR        dm6/dm6_reverse 	        Shortcuts/dm6/dm6.fasta	                      >>dm6.buildIndex.runLog	2>&1



mkdir  -p  OtherGenomes/Ecoli                					
lastdb -uBISF        OtherGenomes/Ecoli/Ecoli_forward        Shortcuts/OtherGenomes/Ecoli.fasta	      >>Ecoli.buildIndex.runLog	2>&1
lastdb -uBISR        OtherGenomes/Ecoli/Ecoli_reverse        Shortcuts/OtherGenomes/Ecoli.fasta	      >>Ecoli.buildIndex.runLog	2>&1

mkdir -p OtherGenomes/phiX174
lastdb -uBISF        OtherGenomes/phiX174/phiX174_forward	        Shortcuts/OtherGenomes/phiX174.fasta	                      >>phiX174.buildIndex.runLog	2>&1
lastdb -uBISR        OtherGenomes/phiX174/phiX174_reverse	        Shortcuts/OtherGenomes/phiX174.fasta	                      >>phiX174.buildIndex.runLog	2>&1

mkdir -p OtherGenomes/Adapters
lastdb -uBISF        OtherGenomes/Adapters/Adapters_forward 	        Shortcuts/OtherGenomes/Adapters.fasta	                      >>Adapters.buildIndex.runLog	2>&1
lastdb -uBISR        OtherGenomes/Adapters/Adapters_reverse 	        Shortcuts/OtherGenomes/Adapters.fasta	                      >>Adapters.buildIndex.runLog	2>&1


mkdir -p OtherGenomes/UniVectors
lastdb -uBISF        OtherGenomes/UniVectors/UniVectors_forward 	        Shortcuts/OtherGenomes/UniVectors.fasta	                      >>UniVectors.buildIndex.runLog	2>&1
lastdb -uBISR        OtherGenomes/UniVectors/UniVectors_reverse 	        Shortcuts/OtherGenomes/UniVectors.fasta	                      >>UniVectors.buildIndex.runLog	2>&1




