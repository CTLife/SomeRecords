out=8-merge-motif/1-bed
mkdir -p  $out


DRACH_TP_1=4-DRACH-TP/5-motif/1-bed/3UTR.bed
DRACH_TP_2=4-DRACH-TP/5-motif/1-bed/5UTR.bed
DRACH_TP_3=4-DRACH-TP/5-motif/1-bed/all.bed
DRACH_TP_4=4-DRACH-TP/5-motif/1-bed/Exon.bed

nonDRACH_TP_1=6-nonDRACH-TP/5-motif/1-bed/3UTR.bed
nonDRACH_TP_2=6-nonDRACH-TP/5-motif/1-bed/5UTR.bed
nonDRACH_TP_3=6-nonDRACH-TP/5-motif/1-bed/all.bed
nonDRACH_TP_4=6-nonDRACH-TP/5-motif/1-bed/Exon.bed
    
cat $DRACH_TP_1   $nonDRACH_TP_1  > $out/TP.3UTR.bed
cat $DRACH_TP_2   $nonDRACH_TP_2  > $out/TP.5UTR.bed
cat $DRACH_TP_3   $nonDRACH_TP_3  > $out/TP.all.bed
cat $DRACH_TP_4   $nonDRACH_TP_4  > $out/TP.Exon.bed




DRACH_FP_1=5-DRACH-FP/5-motif/1-bed/3UTR.bed
DRACH_FP_2=5-DRACH-FP/5-motif/1-bed/5UTR.bed
DRACH_FP_3=5-DRACH-FP/5-motif/1-bed/all.bed
DRACH_FP_4=5-DRACH-FP/5-motif/1-bed/Exon.bed

nonDRACH_FP_1=7-nonDRACH-FP/5-motif/1-bed/3UTR.bed
nonDRACH_FP_2=7-nonDRACH-FP/5-motif/1-bed/5UTR.bed
nonDRACH_FP_3=7-nonDRACH-FP/5-motif/1-bed/all.bed
nonDRACH_FP_4=7-nonDRACH-FP/5-motif/1-bed/Exon.bed
    
cat $DRACH_FP_1   $nonDRACH_FP_1  > $out/FP.3UTR.bed
cat $DRACH_FP_2   $nonDRACH_FP_2  > $out/FP.5UTR.bed
cat $DRACH_FP_3   $nonDRACH_FP_3  > $out/FP.all.bed
cat $DRACH_FP_4   $nonDRACH_FP_4  > $out/FP.Exon.bed






