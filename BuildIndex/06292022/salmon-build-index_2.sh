
mygenome="../../../14_Genomes/Gencode/Mouse/GRCm39.primary_assembly.genome.fa"
mytransc="../../../14_Genomes/Gencode/Mouse/gencode.vM29.transcripts.fa" 
mydir="mm39"
mkdir  $mydir
grep "^>"   $mygenome   | cut -d  " " -f 1 > $mydir/decoys.txt
sed -i.bak -e 's/>//g' $mydir/decoys.txt
cat  $mytransc   $mygenome  > $mydir/gentrome.fa
salmon index   -t $mydir/gentrome.fa   -d $mydir/decoys.txt   -p 12    -i $mydir/salmon_index   --gencode   -k 31
 

 


