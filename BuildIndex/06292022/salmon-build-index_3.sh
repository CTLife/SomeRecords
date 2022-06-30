
mygenome="../../../14_Genomes/UCSC/Green_monkey/chlSab2.fa"
mytransc="../../../14_Genomes/UCSC/Green_monkey/mrna.fa" 
mydir="chlSab2"
mkdir  $mydir
grep "^>"   $mygenome   | cut -d  " " -f 1 > $mydir/decoys.txt
sed -i.bak -e 's/>//g' $mydir/decoys.txt
cat  $mytransc   $mygenome  > $mydir/gentrome.fa
salmon index   -t $mydir/gentrome.fa   -d $mydir/decoys.txt   -p 12    -i $mydir/salmon_index   -k 31
 

 


mygenome="../../../14_Genomes/UCSC/Mouse_lemur/micMur2.fa"
mytransc="../../../14_Genomes/UCSC/Mouse_lemur/mrna.fa" 
mydir="micMur2"
mkdir  $mydir
grep "^>"   $mygenome   | cut -d  " " -f 1 > $mydir/decoys.txt
sed -i.bak -e 's/>//g' $mydir/decoys.txt
cat  $mytransc   $mygenome  > $mydir/gentrome.fa
salmon index   -t $mydir/gentrome.fa   -d $mydir/decoys.txt   -p 12    -i $mydir/salmon_index   -k 31
 

 


mygenome="../../../14_Genomes/UCSC/Zebrafish/danRer11.fa"
mytransc="../../../14_Genomes/UCSC/Zebrafish/mrna.fa" 
mydir="danRer11"
mkdir  $mydir
grep "^>"   $mygenome   | cut -d  " " -f 1 > $mydir/decoys.txt
sed -i.bak -e 's/>//g' $mydir/decoys.txt
cat  $mytransc   $mygenome  > $mydir/gentrome.fa
salmon index   -t $mydir/gentrome.fa   -d $mydir/decoys.txt   -p 12    -i $mydir/salmon_index   -k 31
 

 


