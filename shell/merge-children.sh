outpath="merge-children"
outpath2="$outpath/peaks"
mkdir -p $outpath
mkdir -p $outpath2
mergePeaks  -d  given  -prefix $outpath2/H3K4me1  -matrix $outpath/matrix   -venn $outpath/venn.txt    ART*  NC*   > $outpath/runLog.txt 2>&1  


            
