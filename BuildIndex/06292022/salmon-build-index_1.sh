
mygenome="../../../14_Genomes/Gencode/Human/GRCh38.primary_assembly.genome.fa"
mytransc="../../../14_Genomes/Gencode/Human/gencode.v40.transcripts.fa"
mydir="hg38"
mkdir  $mydir
grep "^>"   $mygenome   | cut -d  " " -f 1 > $mydir/decoys.txt
sed -i.bak -e 's/>//g' $mydir/decoys.txt
cat  $mytransc   $mygenome  > $mydir/gentrome.fa
salmon index   -t $mydir/gentrome.fa   -d $mydir/decoys.txt   -p 12    -i $mydir/salmon_index   --gencode   -k 31



