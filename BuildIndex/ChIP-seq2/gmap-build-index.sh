mkdir    sacCer3
gmap_build   --dir=sacCer3	 --db=sacCer3       Shortcuts/sacCer3/sacCer3.fasta      >> sacCer3.runLog	 2>&1          

mkdir    hg38
gmap_build   --dir=hg38	 --db=hg38       Shortcuts/hg38/hg38.fasta      >> hg38.runLog	 2>&1          

mkdir    mm10
gmap_build   --dir=mm10	 --db=mm10       Shortcuts/mm10/mm10.fasta      >> mm10.runLog	 2>&1          



mkdir    ce11
gmap_build   --dir=ce11	 --db=ce11       Shortcuts/ce11/ce11.fasta      >> ce11.runLog	 2>&1          

mkdir    danRer10
gmap_build   --dir=danRer10	 --db=danRer10       Shortcuts/danRer10/danRer10.fasta      >> danRer10.runLog	 2>&1          

mkdir    dm6
gmap_build   --dir=dm6	 --db=dm6       Shortcuts/dm6/dm6.fasta      >> dm6.runLog	 2>&1          





mkdir    OtherGenomes/Adapters
gmap_build   --dir=OtherGenomes/Adapters	 --db=Adapters       Shortcuts/OtherGenomes/Adapters.fasta      >> Adapters.runLog	 2>&1          

mkdir    OtherGenomes/Ecoli
gmap_build   --dir=OtherGenomes/Ecoli	 --db=Ecoli       Shortcuts/OtherGenomes/Ecoli.fasta      >> Ecoli.runLog	 2>&1          

mkdir    OtherGenomes/phiX174
gmap_build   --dir=OtherGenomes/phiX174	 --db=phiX174       Shortcuts/OtherGenomes/phiX174.fasta      >> phiX174.runLog	 2>&1          

 
mkdir    OtherGenomes/UniVectors
gmap_build   --dir=OtherGenomes/UniVectors	 --db=UniVectors       Shortcuts/OtherGenomes/UniVectors.fasta      >> UniVectors.runLog	 2>&1          

 





mkdir    hg38.ensembl
gmap_build   --dir=hg38.ensembl	 --db=hg38.ensembl       Shortcuts/hg38/hg38.ensembl.genome.fasta      >> hg38.ensembl.runLog	 2>&1          

mkdir    mm10.ensembl
gmap_build   --dir=mm10.ensembl	 --db=mm10.ensembl       Shortcuts/mm10/mm10.ensembl.genome.fasta      >> mm10.ensembl.runLog	 2>&1          












