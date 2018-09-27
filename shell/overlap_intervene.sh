file1_hyper="1_NC_vs_IVF-fresh/twins_hyper.bed"
file1_hypo="1_NC_vs_IVF-fresh/twins_hypo.bed"
file2_hyper="2_NC_vs_ICSI-fresh/twins_hyper.bed"
file2_hypo="2_NC_vs_ICSI-fresh/twins_hypo.bed"
file3_hyper="3_NC_vs_IVF-frozen/twins_hyper.bed"
file3_hypo="3_NC_vs_IVF-frozen/twins_hypo.bed"
file4_hyper="4_NC_vs_ICSI-frozen/twins_hyper.bed"
file4_hypo="4_NC_vs_ICSI-frozen/twins_hypo.bed"

outPath1="100_4groups_overlap/hyper"
outPath2="100_4groups_overlap/hypo"
mkdir -p $outPath1
mkdir -p $outPath2

intervene  venn  --input   $file1_hyper   $file2_hyper   $file3_hyper   $file4_hyper  --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_ICSI-fresh,NC_vs_IVF-frozen,NC_vs_ICSI-frozen             \
--output $outPath1          --save-overlaps      --dpi 1200    --figtype svg     
    
intervene  venn  --input   $file1_hypo   $file2_hypo   $file3_hypo   $file4_hypo  --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_ICSI-fresh,NC_vs_IVF-frozen,NC_vs_ICSI-frozen             \
--output $outPath2          --save-overlaps      --dpi 1200    --figtype svg     




###############################################################################
outPath3_1_hyper="101_2groups_overlap/1_hyper"
outPath3_1_hypo="101_2groups_overlap/1_hypo"
mkdir -p $outPath3_1_hyper
mkdir -p $outPath3_1_hypo
intervene  venn  --input   $file1_hyper   $file2_hyper     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_ICSI-fresh              \
--output $outPath3_1_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
intervene  venn  --input   $file1_hypo   $file2_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_ICSI-fresh              \
--output $outPath3_1_hypo      --save-overlaps      --dpi 1200    --figtype svg     
   




outPath3_2_hyper="101_2groups_overlap/2_hyper"
outPath3_2_hypo="101_2groups_overlap/2_hypo"
mkdir -p $outPath3_2_hyper
mkdir -p $outPath3_2_hypo
intervene  venn  --input   $file1_hyper   $file3_hyper     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_IVF-frozen             \
--output $outPath3_2_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
intervene  venn  --input   $file1_hypo   $file3_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_IVF-frozen             \
--output $outPath3_2_hypo      --save-overlaps      --dpi 1200    --figtype svg     





outPath3_3_hyper="101_2groups_overlap/3_hyper"
outPath3_3_hypo="101_2groups_overlap/3_hypo"
mkdir -p $outPath3_3_hyper
mkdir -p $outPath3_3_hypo
intervene  venn  --input   $file1_hyper   $file4_hyper     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_ICSI-frozen             \
--output $outPath3_3_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
intervene  venn  --input   $file1_hypo   $file4_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh,NC_vs_ICSI-frozen             \
--output $outPath3_3_hypo      --save-overlaps      --dpi 1200    --figtype svg     
     




outPath3_4_hyper="101_2groups_overlap/4_hyper"
outPath3_4_hypo="101_2groups_overlap/4_hypo"
mkdir -p $outPath3_4_hyper
mkdir -p $outPath3_4_hypo
intervene  venn  --input   $file2_hyper   $file3_hyper     --bedtools-options f=0.5,e      \
--names=NC_vs_ICSI-fresh,NC_vs_IVF-frozen             \
--output $outPath3_4_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
intervene  venn  --input   $file2_hypo   $file3_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_ICSI-fresh,NC_vs_IVF-frozen             \
--output $outPath3_4_hypo      --save-overlaps      --dpi 1200    --figtype svg     
 
     




outPath3_5_hyper="101_2groups_overlap/5_hyper"
outPath3_5_hypo="101_2groups_overlap/5_hypo"
mkdir -p $outPath3_5_hyper
mkdir -p $outPath3_5_hypo
intervene  venn  --input   $file2_hyper   $file4_hyper     --bedtools-options f=0.5,e      \
--names=NC_vs_ICSI-fresh,NC_vs_ICSI-frozen             \
--output $outPath3_5_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
intervene  venn  --input   $file2_hypo   $file4_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_ICSI-fresh,NC_vs_ICSI-frozen             \
--output $outPath3_5_hypo      --save-overlaps      --dpi 1200    --figtype svg     
 
 
     




outPath3_6_hyper="101_2groups_overlap/6_hyper"
outPath3_6_hypo="101_2groups_overlap/6_hypo"
mkdir -p $outPath3_6_hyper
mkdir -p $outPath3_6_hypo
intervene  venn  --input   $file3_hyper   $file4_hyper     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-frozen,NC_vs_ICSI-frozen             \
--output $outPath3_6_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
intervene  venn  --input   $file3_hypo   $file4_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-frozen,NC_vs_ICSI-frozen             \
--output $outPath3_6_hypo      --save-overlaps      --dpi 1200    --figtype svg     
    
   




 



###############################################################################
outPath4_1_hyper="102_4groups_overlap/1"
mkdir -p $outPath4_1_hyper
intervene  venn  --input   $file1_hyper   $file2_hyper     $file1_hypo   $file2_hypo   --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh_hyper,NC_vs_ICSI-fresh_hyper,NC_vs_IVF-fresh_hypo,NC_vs_ICSI-fresh_hypo               \
--output $outPath4_1_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
  




outPath4_2_hyper="102_4groups_overlap/2"
mkdir -p $outPath4_2_hyper
intervene  venn  --input   $file1_hyper   $file3_hyper    $file1_hypo   $file3_hypo      --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh_hyper,NC_vs_IVF-frozen_hyper,NC_vs_IVF-fresh_hypo,NC_vs_IVF-frozen_hypo             \
--output $outPath4_2_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
 




outPath4_3_hyper="102_4groups_overlap/3"
mkdir -p $outPath4_3_hyper
intervene  venn  --input   $file1_hyper   $file4_hyper   $file1_hypo   $file4_hypo    --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-fresh_hyper,NC_vs_ICSI-frozen_hyper,NC_vs_IVF-fresh_hypo,NC_vs_ICSI-frozen_hypo             \
--output $outPath4_3_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
 




outPath4_4_hyper="102_4groups_overlap/4"
mkdir -p $outPath4_4_hyper
intervene  venn  --input   $file2_hyper   $file3_hyper    $file2_hypo   $file3_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_ICSI-fresh_hyper,NC_vs_IVF-frozen_hyper,NC_vs_ICSI-fresh_hypo,NC_vs_IVF-frozen_hypo             \
--output $outPath4_4_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
 
     




outPath4_5_hyper="102_4groups_overlap/5"
mkdir -p $outPath4_5_hyper
intervene  venn  --input   $file2_hyper   $file4_hyper   $file2_hypo   $file4_hypo     --bedtools-options f=0.5,e      \
--names=NC_vs_ICSI-fresh_hyper,NC_vs_ICSI-frozen_hyper,NC_vs_ICSI-fresh_hypo,NC_vs_ICSI-frozen_hypo             \
--output $outPath4_5_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
 
 
     




outPath4_6_hyper="102_4groups_overlap/6"
mkdir -p $outPath4_6_hyper
intervene  venn  --input   $file3_hyper   $file4_hyper    $file3_hypo   $file4_hypo    --bedtools-options f=0.5,e      \
--names=NC_vs_IVF-frozen_hyper,NC_vs_ICSI-frozen_hyper,NC_vs_IVF-frozen_hypo,NC_vs_ICSI-frozen_hypo             \
--output $outPath4_6_hyper      --save-overlaps      --dpi 1200    --figtype svg     
   
    
   




 
