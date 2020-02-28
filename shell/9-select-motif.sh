indir=8-merge-motif/1-bed
out=9-select-motif/1-bed
mkdir -p  $out

cat  $indir/TP.3UTR.bed   | awk  -F "\t"  '$5>=10{print $0}'  > $out/TP.3UTR.bed
cat  $indir/TP.5UTR.bed   | awk  -F "\t"  '$5>=10{print $0}'  > $out/TP.5UTR.bed
cat  $indir/TP.all.bed    | awk  -F "\t"  '$5>=10{print $0}'  > $out/TP.all.bed
cat  $indir/TP.Exon.bed   | awk  -F "\t"  '$5>=10{print $0}'  > $out/TP.Exon.bed
   
cat  $indir/FP.3UTR.bed   | awk  -F "\t"  '$5>=10{print $0}'  > $out/FP.3UTR.bed
cat  $indir/FP.5UTR.bed   | awk  -F "\t"  '$5>=10{print $0}'  > $out/FP.5UTR.bed
cat  $indir/FP.all.bed    | awk  -F "\t"  '$5>=10{print $0}'  > $out/FP.all.bed
cat  $indir/FP.Exon.bed   | awk  -F "\t"  '$5>=10{print $0}'  > $out/FP.Exon.bed






