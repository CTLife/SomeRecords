#    venn                Venn diagram of intersection of genomic regions or list sets (upto 6-way).
#    upset               UpSet diagram of intersection of genomic regions or list sets.
#    pairwise            Pairwise intersection and heatmap of N genomic region sets in <BED/GTF/GFF> format.




A0=2-finalBED/1_BA9.bed
A1=2-finalBED/2_BA24.bed
A2=2-finalBED/3_C.bed
A3=2-finalBED/4_H.bed
A4=2-finalBED/5_T.bed




outdir="5-intervenn-all"
mkdir  $outdir 






## venn
##############################################################################
outpath1=$outdir/venn
mkdir  $outpath1
intervene  venn  --input  \
$A0  \
$A1  \
$A2  \
$A3  \
$A4  \
--type genomic    --save-overlaps   --output $outpath1  \
--dpi 1200  --figtype pdf   --figsize 12 12  --fontsize 10
##############################################################################




## upset
##############################################################################
outpath2=$outdir/upset
mkdir  $outpath2
intervene  upset  --input  \
$A0  \
$A1  \
$A2  \
$A3  \
$A4  \
--type   genomic    --save-overlaps   --output $outpath2  --ninter 100       \
--dpi 1200  --figtype pdf   --figsize 12 6   
##############################################################################






## pairwise
##############################################################################
outpath1="$outdir/pairwise.fraction"
mkdir -p $outpath1

intervene  pairwise  \
--input  \
$A0  \
$A1  \
$A2  \
$A3  \
$A4  \
--compute  frac  --htype color   --output $outpath1  \
--dpi 1200  --figtype svg   --figsize 10 10  --fontsize 10



outpath2="$outdir/pairwise.count"
mkdir -p $outpath2

intervene  pairwise  \
--input  \
$A0  \
$A1  \
$A2  \
$A3  \
$A4  \
--compute  count  --htype color  --output $outpath2  \
--dpi 1200  --figtype svg   --figsize 10 10   --fontsize 10
##############################################################################


























outpath1="3-overlap/pairwise.fraction"
mkdir -p $outpath1

intervene  pairwise  \
--input  \
$A0  \
$A1  \
$A2  \
$A3  \
$A4  \
--compute  frac  --htype color   --output $outpath1 --save-overlaps   \
--ninter  100  --dpi 1200  --figtype svg   --figsize 10 10  --fontsize 10









outpath2="3-overlap/pairwise.count"
mkdir -p $outpath2

intervene  pairwise  \
--input  \
$A0  \
$A1  \
$A2  \
$A3  \
$A4  \
--compute  count  --htype color  --output $outpath2  \
--dpi 1200  --figtype svg   --figsize 10 10   --fontsize 10







