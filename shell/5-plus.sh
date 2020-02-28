out=5-motif/1-bed
mkdir -p  $out


cat 2-select_5percent/3UTR.txt  | awk  '{print $1"\t"$2-6"\t"$3+5"\t"$6"\t"$7"\t"$8}' > $out/3UTR.bed            
cat 2-select_5percent/5UTR.txt  | awk  '{print $1"\t"$2-6"\t"$3+5"\t"$6"\t"$7"\t"$8}' > $out/5UTR.bed            
cat 2-select_5percent/Exon.txt  | awk  '{print $1"\t"$2-6"\t"$3+5"\t"$6"\t"$7"\t"$8}' > $out/Exon.bed            
cat $out/*.bed > $out/all.bed


