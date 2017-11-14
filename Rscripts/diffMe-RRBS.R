###########################################################################
## 1. Read the raw coverage files.
###########################################################################
myOutDir <- "100-diffMe-1"

myFileLists <- list(
"1-Coverage-CpG/61_NC-BS2-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/61_NC-BS2-D-Girl-merge_Rep1.bismark.cov",

"1-Coverage-CpG/66_ART-E29-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/66_ART-E29-D-Girl-merge_Rep1.bismark.cov"
)

mySampleID <- list( "NC1", "NC2",  "IVF-fresh1", "IVF-fresh2"  )
myTreatment <- c( 0, 0, 1, 1)       ## This option will determine the result of unite.

myOutDir_sub1 = paste(myOutDir, "/1-ReadRawFiles",  sep="");
if( ! file.exists(myOutDir_sub1) ) { dir.create(myOutDir_sub1, recursive = TRUE) }

print( myTreatment )
print( length(myFileLists) )
print( length(mySampleID) )
print( length(myTreatment) )


library(methylKit)
library(genomation)
library(ggplot2) 
library(ggfortify)
library(cluster)
library(lfda)
library(MASS)
library(factoextra)
library(magrittr)  
library(dplyr)  
library(rgl)
library(ggbiplot)
library(gdata)

continue_on_error <- function()  {
      print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error) 


MyTheme_1 <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all text elements.           "serif" for a serif font. 所有文本相关属性.
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    ## aspect.ratio = 1,   ##高宽比
    
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis label (element_text; inherits from axis.title)
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis label (element_text; inherits from axis.title)
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 

MySaveGgplot2_1 <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps)         
}



# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1, "2-theLog-of-read-AllTheFiles.txt", sep="/") )
myobj=methRead(myFileLists,
               sample.id=mySampleID,
               assembly="hg38",
               treatment=myTreatment,
               context="CpG",
               pipeline = "bismarkCoverage",
               mincov = 5,       ## >= n
               header = FALSE
)
sink()


sink( file=paste(myOutDir_sub1, "3-all-rawFiles.txt", sep="/") )
    print(myFileLists)
    print("#########################")
    print("#########################")
    print(myobj)
sink()


sink( file=paste(myOutDir_sub1, "4A-all-files-1reads.txt", sep="/") )
    print(myobj)
sink()




#Filtering samples based on read coverage
filtered.myobj_1 = filterByCoverage(myobj,  lo.count=5,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)
sink( file=paste(myOutDir_sub1, "4B-all-files-5reads.txt", sep="/") )
    print(filtered.myobj_1)
sink()

filtered.myobj_2 = filterByCoverage(myobj,  lo.count=5,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.9)
sink( file=paste(myOutDir_sub1, "4C-all-files-5reads-rmHigh.txt", sep="/") )
print(filtered.myobj_2)
sink()

filtered.myobj_3 = filterByCoverage(myobj,  lo.count=10,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)
sink( file=paste(myOutDir_sub1, "4D-all-files-10reads.txt", sep="/") )
print(filtered.myobj_3)
sink()

filtered.myobj_4 = filterByCoverage(myobj,  lo.count=10,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.9)
sink( file=paste(myOutDir_sub1, "4E-all-files-10reads-rmHigh.txt", sep="/") )
print(filtered.myobj_4)
sink()





sink( file=paste(myOutDir_sub1, "5-dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print( "######################" )
  print(   myFileLists[[i]]  )
  print(   dim(myobj[[i]])  )
  print(   dim(filtered.myobj_1[[i]])  )
  print(   dim(filtered.myobj_2[[i]])  )
  print(   dim(filtered.myobj_3[[i]])  )
  print(   dim(filtered.myobj_4[[i]])  )
}
sink()



## Merging samples
## Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. 
## This provides better coverage, but only advised when looking at CpG  
## methylation (for CpH methylation this will cause wrong results in subsequent analyses).

sink(file=paste(myOutDir_sub1, "6-log-merged-overlapSites.txt", sep="/") )

meth_1 = unite(filtered.myobj_1, destrand=FALSE  )
head(meth_1)
dim(meth_1)

meth_2 = unite(filtered.myobj_2, destrand=FALSE  )
head(meth_2)
dim(meth_2)

meth_3 = unite(filtered.myobj_3, destrand=FALSE  )
head(meth_3)
dim(meth_3)

meth_4 = unite(filtered.myobj_4, destrand=FALSE  )
head(meth_4)
dim(meth_4)

sink()




sink(file=paste(myOutDir_sub1, "7-dimensions-merged-overlap.txt", sep="/") )
    print( dim(meth_1) )
    print( dim(meth_2) )
    print( dim(meth_3) )
    print( dim(meth_4) )
sink()





sink(file=paste(myOutDir_sub1, "8-log-percMethylation.txt", sep="/") )

mat_1 = percMethylation(meth_1)
head(mat_1)
dim(mat_1)

mat_2 = percMethylation(meth_2)
head(mat_2)
dim(mat_2)

mat_3 = percMethylation(meth_3)
head(mat_3)
dim(mat_3)

mat_4 = percMethylation(meth_4)
head(mat_4)
dim(mat_4)

sink()



sink(file=paste(myOutDir_sub1, "9-dimensions-methylationMatrix.txt", sep="/") )
  print( dim(mat_1) )
  print( dim(mat_2) )
  print( dim(mat_3) )
  print( dim(mat_4) )
sink()




write.table(mat_1 , 
            file = paste(myOutDir_sub1,"10B_mathylationLevel_1_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(mat_2 , 
            file = paste(myOutDir_sub1,"10C_mathylationLevel_2_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3 , 
            file = paste(myOutDir_sub1,"10D_mathylationLevel_3_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4 , 
            file = paste(myOutDir_sub1,"10E_mathylationLevel_4_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

###########################################################################
###########################################################################











######################################################################################################################################################
######################################################################################################################################################
myOutDir_sub3 = paste(myOutDir, "/3-Cov-5reads-noRemoveHigh",  sep="");
if( ! file.exists(myOutDir_sub3) ) { dir.create(myOutDir_sub3, recursive = TRUE) }

pdf( file=paste(myOutDir_sub3, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_1[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub3, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_1[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub3, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_1[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub3, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_1[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub3, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
    clusterSamples(meth_1, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_1, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_1, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_1, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_sub3, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1, screeplot=TRUE)
PCASamples(meth_1)
dev.off()




##################
myDiff_1 = calculateDiffMeth(meth_1, num.cores=8)

myDiff25p.hypo_1  = getMethylDiff(myDiff_1, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p.hyper_1 = getMethylDiff(myDiff_1, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p_1 = getMethylDiff(myDiff_1, difference=25, qvalue=0.05)


write.table(myDiff_1 , file = paste(myOutDir_sub3,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1 , file = paste(myOutDir_sub3,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1 , file = paste(myOutDir_sub3,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1 , file = paste(myOutDir_sub3,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub3, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub3, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_1,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions
myRefSeqGenes = "/home/yp/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands  = "/home/yp/AnnotationBED/hg38_CpG_islands.bed"
myRepeats     = "/home/yp/AnnotationBED/hg38_Repeats_rmsk.bed"


diffGeneAnn_1_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub3, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1,"GRanges"),
                                    cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                    feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub3, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1,"GRanges"),
                                             repeat.obj_1$Repeats,  repeat.obj_1$others,
                                             feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub3, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub3, "8A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub3, "8B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1,"GRanges"),
                                              cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                              feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub3, "8C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub3, "8D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()





diffrepeatann_1_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1,"GRanges"),
                                                repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub3, "8E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "8F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()












### Tiling windows analysis
tiles_1=tileMethylCounts(filtered.myobj_1, win.size=1000, step.size=1000)
head(tiles_1[[1]], 3)

meth_1_region = unite(tiles_1, destrand=FALSE  )
head(meth_1_region)
dim(meth_1_region)

mat_1_region = percMethylation(meth_1_region)
head(mat_1_region)
dim(mat_1_region)

write.table(meth_1_region , 
            file = paste(myOutDir_sub3,"10A_meth_1_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1_region , 
            file = paste(myOutDir_sub3,"10B_mat_1_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_1_region = calculateDiffMeth(meth_1_region, num.cores=8)

myDiff25p.hypo_1_region  = getMethylDiff(myDiff_1_region, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p.hyper_1_region = getMethylDiff(myDiff_1_region, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p_1_region       = getMethylDiff(myDiff_1_region, difference=25, qvalue=0.05)


write.table(myDiff_1_region , file = paste(myOutDir_sub3,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1_region , file = paste(myOutDir_sub3,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1_region , file = paste(myOutDir_sub3,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1_region , file = paste(myOutDir_sub3,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub3, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1_region,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub3, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_1_region,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions

diffGeneAnn_1_hypo_region = annotateWithGeneParts(as(myDiff25p.hypo_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub3, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo_region)
sink()
pdf( file=paste(myOutDir_sub3, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_1_region,"GRanges"),
                                             cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                             feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub3, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo_region,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_1_region,"GRanges"),
                                                repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub3, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper_region = annotateWithGeneParts(as(myDiff25p.hyper_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub3, "14A-distribution-onGenes-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper_region,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub3, "14B-distribution-onGenes-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p.hyper_1_region,"GRanges"),
                                              cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                              feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub3, "14C-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper_region,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub3, "14D-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


diffrepeatann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p.hyper_1_region,"GRanges"),
                                                 repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                 feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub3, "14E-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub3, "14F-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


















######################################################################################################################################################
######################################################################################################################################################
myOutDir_sub4 = paste(myOutDir, "/4-Cov-5reads-RemoveHigh",  sep="");
if( ! file.exists(myOutDir_sub4) ) { dir.create(myOutDir_sub4, recursive = TRUE) }

pdf( file=paste(myOutDir_sub4, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_2[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub4, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_2[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub4, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_2[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub4, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_2[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub4, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_sub4, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2, screeplot=TRUE)
PCASamples(meth_2)
dev.off()




##################
myDiff_2 = calculateDiffMeth(meth_2, num.cores=8)

myDiff25p2.hypo_1  = getMethylDiff(myDiff_2, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p2.hyper_1 = getMethylDiff(myDiff_2, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p2_1 = getMethylDiff(myDiff_2, difference=25, qvalue=0.05)


write.table(myDiff_2 , file = paste(myOutDir_sub4,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p2.hypo_1 , file = paste(myOutDir_sub4,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p2.hyper_1 , file = paste(myOutDir_sub4,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p2_1 , file = paste(myOutDir_sub4,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub4, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub4, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_2,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions
myRefSeqGenes = "/home/yp/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands  = "/home/yp/AnnotationBED/hg38_CpG_islands.bed"
myRepeats     = "/home/yp/AnnotationBED/hg38_Repeats_rmsk.bed"

gene.obj_1 = readTranscriptFeatures(myRefSeqGenes)
cpg.obj_1  = readFeatureFlank(myCpGIslands, feature.flank.name=c("CpGi","shores"))
repeat.obj_1  = readFeatureFlank(myRepeats, feature.flank.name=c("Repeats","others"))

diffGeneAnn_1_hypo = annotateWithGeneParts(as(myDiff25p2.hypo_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub4, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo = annotateWithFeatureFlank(as(myDiff25p2.hypo_1,"GRanges"),
                                             cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                             feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub4, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo = annotateWithFeatureFlank(as(myDiff25p2.hypo_1,"GRanges"),
                                                repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub4, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper = annotateWithGeneParts(as(myDiff25p2.hyper_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub4, "8A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub4, "8B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper = annotateWithFeatureFlank(as(myDiff25p2.hyper_1,"GRanges"),
                                              cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                              feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub4, "8C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub4, "8D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()





diffrepeatann_1_hyper = annotateWithFeatureFlank(as(myDiff25p2.hyper_1,"GRanges"),
                                                 repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                 feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub4, "8E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "8F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()












### Tiling windows analysis
tiles_1=tileMethylCounts(filtered.myobj_2, win.size=1000, step.size=1000)
head(tiles_1[[1]], 3)

meth_2_region = unite(tiles_1, destrand=FALSE  )
head(meth_2_region)
dim(meth_2_region)

mat_1_region = percMethylation(meth_2_region)
head(mat_1_region)
dim(mat_1_region)

write.table(meth_2_region , 
            file = paste(myOutDir_sub4,"10A_meth_2_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1_region , 
            file = paste(myOutDir_sub4,"10B_mat_1_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_2_region = calculateDiffMeth(meth_2_region, num.cores=8)

myDiff25p2.hypo_1_region  = getMethylDiff(myDiff_2_region, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p2.hyper_1_region = getMethylDiff(myDiff_2_region, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p2_1_region       = getMethylDiff(myDiff_2_region, difference=25, qvalue=0.05)


write.table(myDiff_2_region , file = paste(myOutDir_sub4,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p2.hypo_1_region , file = paste(myOutDir_sub4,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p2.hyper_1_region , file = paste(myOutDir_sub4,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p2_1_region , file = paste(myOutDir_sub4,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub4, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2_region,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub4, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_2_region,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions

diffGeneAnn_1_hypo_region = annotateWithGeneParts(as(myDiff25p2.hypo_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub4, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo_region)
sink()
pdf( file=paste(myOutDir_sub4, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p2.hypo_1_region,"GRanges"),
                                                    cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                                    feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub4, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo_region,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p2.hypo_1_region,"GRanges"),
                                                       repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                       feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub4, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper_region = annotateWithGeneParts(as(myDiff25p2.hyper_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub4, "14A-distribution-onGenes-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper_region,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub4, "14B-distribution-onGenes-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p2.hyper_1_region,"GRanges"),
                                                     cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                                     feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub4, "14C-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper_region,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub4, "14D-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


diffrepeatann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p2.hyper_1_region,"GRanges"),
                                                        repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                        feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub4, "14E-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub4, "14F-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()














######################################################################################################################################################
######################################################################################################################################################
myOutDir_sub5 = paste(myOutDir, "/5-Cov-10reads-noRemoveHigh",  sep="");
if( ! file.exists(myOutDir_sub5) ) { dir.create(myOutDir_sub5, recursive = TRUE) }

pdf( file=paste(myOutDir_sub5, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_3[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub5, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_3[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub5, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_3[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub5, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_3[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub5, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_sub5, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3, screeplot=TRUE)
PCASamples(meth_3)
dev.off()




##################
myDiff_3 = calculateDiffMeth(meth_3, num.cores=8)

myDiff25p3.hypo_1  = getMethylDiff(myDiff_3, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p3.hyper_1 = getMethylDiff(myDiff_3, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p3_1 = getMethylDiff(myDiff_3, difference=25, qvalue=0.05)


write.table(myDiff_3 , file = paste(myOutDir_sub5,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p3.hypo_1 , file = paste(myOutDir_sub5,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p3.hyper_1 , file = paste(myOutDir_sub5,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p3_1 , file = paste(myOutDir_sub5,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub5, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub5, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_3,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions
myRefSeqGenes = "/home/yp/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands  = "/home/yp/AnnotationBED/hg38_CpG_islands.bed"
myRepeats     = "/home/yp/AnnotationBED/hg38_Repeats_rmsk.bed"

 
diffGeneAnn_1_hypo = annotateWithGeneParts(as(myDiff25p3.hypo_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub5, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo = annotateWithFeatureFlank(as(myDiff25p3.hypo_1,"GRanges"),
                                             cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                             feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub5, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo = annotateWithFeatureFlank(as(myDiff25p3.hypo_1,"GRanges"),
                                                repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub5, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper = annotateWithGeneParts(as(myDiff25p3.hyper_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub5, "8A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub5, "8B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper = annotateWithFeatureFlank(as(myDiff25p3.hyper_1,"GRanges"),
                                              cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                              feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub5, "8C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub5, "8D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()





diffrepeatann_1_hyper = annotateWithFeatureFlank(as(myDiff25p3.hyper_1,"GRanges"),
                                                 repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                 feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub5, "8E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "8F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()












### Tiling windows analysis
tiles_1=tileMethylCounts(filtered.myobj_3, win.size=1000, step.size=1000)
head(tiles_1[[1]], 3)

meth_3_region = unite(tiles_1, destrand=FALSE  )
head(meth_3_region)
dim(meth_3_region)

mat_1_region = percMethylation(meth_3_region)
head(mat_1_region)
dim(mat_1_region)

write.table(meth_3_region , 
            file = paste(myOutDir_sub5,"10A_meth_3_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1_region , 
            file = paste(myOutDir_sub5,"10B_mat_1_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_3_region = calculateDiffMeth(meth_3_region, num.cores=8)

myDiff25p3.hypo_1_region  = getMethylDiff(myDiff_3_region, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p3.hyper_1_region = getMethylDiff(myDiff_3_region, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p3_1_region       = getMethylDiff(myDiff_3_region, difference=25, qvalue=0.05)


write.table(myDiff_3_region , file = paste(myOutDir_sub5,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p3.hypo_1_region , file = paste(myOutDir_sub5,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p3.hyper_1_region , file = paste(myOutDir_sub5,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p3_1_region , file = paste(myOutDir_sub5,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub5, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3_region,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub5, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_3_region,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions

diffGeneAnn_1_hypo_region = annotateWithGeneParts(as(myDiff25p3.hypo_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub5, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo_region)
sink()
pdf( file=paste(myOutDir_sub5, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p3.hypo_1_region,"GRanges"),
                                                    cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                                    feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub5, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo_region,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p3.hypo_1_region,"GRanges"),
                                                       repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                       feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub5, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper_region = annotateWithGeneParts(as(myDiff25p3.hyper_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub5, "14A-distribution-onGenes-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper_region,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub5, "14B-distribution-onGenes-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p3.hyper_1_region,"GRanges"),
                                                     cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                                     feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub5, "14C-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper_region,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub5, "14D-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


diffrepeatann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p3.hyper_1_region,"GRanges"),
                                                        repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                        feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub5, "14E-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub5, "14F-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


































######################################################################################################################################################
######################################################################################################################################################
myOutDir_sub6 = paste(myOutDir, "/6-Cov-10reads-RemoveHigh",  sep="");
if( ! file.exists(myOutDir_sub6) ) { dir.create(myOutDir_sub6, recursive = TRUE) }

pdf( file=paste(myOutDir_sub6, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_4[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub6, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_4[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub6, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_4[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub6, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_4[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub6, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_sub6, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4, screeplot=TRUE)
PCASamples(meth_4)
dev.off()




##################
myDiff_4 = calculateDiffMeth(meth_4, num.cores=8)

myDiff25p4.hypo_1  = getMethylDiff(myDiff_4, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p4.hyper_1 = getMethylDiff(myDiff_4, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p4_1 = getMethylDiff(myDiff_4, difference=25, qvalue=0.05)


write.table(myDiff_4 , file = paste(myOutDir_sub6,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p4.hypo_1 , file = paste(myOutDir_sub6,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p4.hyper_1 , file = paste(myOutDir_sub6,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p4_1 , file = paste(myOutDir_sub6,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub6, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub6, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_4,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions
myRefSeqGenes = "/home/yp/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands  = "/home/yp/AnnotationBED/hg38_CpG_islands.bed"
myRepeats     = "/home/yp/AnnotationBED/hg38_Repeats_rmsk.bed"

 
diffGeneAnn_1_hypo = annotateWithGeneParts(as(myDiff25p4.hypo_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub6, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo = annotateWithFeatureFlank(as(myDiff25p4.hypo_1,"GRanges"),
                                             cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                             feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub6, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo = annotateWithFeatureFlank(as(myDiff25p4.hypo_1,"GRanges"),
                                                repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub6, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper = annotateWithGeneParts(as(myDiff25p4.hyper_1,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub6, "8A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub6, "8B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper = annotateWithFeatureFlank(as(myDiff25p4.hyper_1,"GRanges"),
                                              cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                              feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub6, "8C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub6, "8D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()





diffrepeatann_1_hyper = annotateWithFeatureFlank(as(myDiff25p4.hyper_1,"GRanges"),
                                                 repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                 feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub6, "8E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "8F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()












### Tiling windows analysis
tiles_1=tileMethylCounts(filtered.myobj_4, win.size=1000, step.size=1000)
head(tiles_1[[1]], 3)

meth_4_region = unite(tiles_1, destrand=FALSE  )
head(meth_4_region)
dim(meth_4_region)

mat_1_region = percMethylation(meth_4_region)
head(mat_1_region)
dim(mat_1_region)

write.table(meth_4_region , 
            file = paste(myOutDir_sub6,"10A_meth_4_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1_region , 
            file = paste(myOutDir_sub6,"10B_mat_1_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_4_region = calculateDiffMeth(meth_4_region, num.cores=8)

myDiff25p4.hypo_1_region  = getMethylDiff(myDiff_4_region, difference=25, qvalue=0.05, type="hypo")  ## less enrich in IVF-fresh
myDiff25p4.hyper_1_region = getMethylDiff(myDiff_4_region, difference=25, qvalue=0.05, type="hyper")  ## more enrich in IVF-fresh
myDiff25p4_1_region       = getMethylDiff(myDiff_4_region, difference=25, qvalue=0.05)


write.table(myDiff_4_region , file = paste(myOutDir_sub6,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p4.hypo_1_region , file = paste(myOutDir_sub6,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p4.hyper_1_region , file = paste(myOutDir_sub6,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p4_1_region , file = paste(myOutDir_sub6,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_sub6, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4_region,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=25)
sink()


pdf( file=paste(myOutDir_sub6, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=10, height=5    )
diffMethPerChr(myDiff_4_region,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=25)
dev.off()


## Annotating differentially methylated bases or regions

diffGeneAnn_1_hypo_region = annotateWithGeneParts(as(myDiff25p4.hypo_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub6, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hypo,  percentage=TRUE)
print(diffGeneAnn_1_hypo_region)
sink()
pdf( file=paste(myOutDir_sub6, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p4.hypo_1_region,"GRanges"),
                                                    cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                                    feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub6, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hypo_region,  percentage=TRUE)
print(diffCpGann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



diffrepeatann_1_hypo_region = annotateWithFeatureFlank(as(myDiff25p4.hypo_1_region,"GRanges"),
                                                       repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                       feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub6, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hypo_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




###############
diffGeneAnn_1_hyper_region = annotateWithGeneParts(as(myDiff25p4.hyper_1_region,"GRanges"),  gene.obj_1)
sink( file=paste(myOutDir_sub6, "14A-distribution-onGenes-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1_hyper_region,  percentage=TRUE)
print(diffGeneAnn_1_hyper)
sink()
pdf( file=paste(myOutDir_sub6, "14B-distribution-onGenes-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1_hyper_region,precedence=TRUE, main="differential methylation annotation")
dev.off()


diffCpGann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p4.hyper_1_region,"GRanges"),
                                                     cpg.obj_1$CpGi,  cpg.obj_1$shores,
                                                     feature.name="CpGi",flank.name="shores")
sink( file=paste(myOutDir_sub6, "14C-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1_hyper_region,  percentage=TRUE)
print(diffCpGann_1_hyper)
sink()
pdf( file=paste(myOutDir_sub6, "14D-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


diffrepeatann_1_hyper_region = annotateWithFeatureFlank(as(myDiff25p4.hyper_1_region,"GRanges"),
                                                        repeat.obj_1$Repeats,  repeat.obj_1$others,
                                                        feature.name="Repeats",flank.name="others")
sink( file=paste(myOutDir_sub6, "14E-distribution-on-hyper_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1_hyper_region,  percentage=TRUE)
print(diffrepeatann_1_hypo)
sink()
pdf( file=paste(myOutDir_sub6, "14F-distribution-onCpGs-hyper_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1_hyper_region, precedence=TRUE, main="differential methylation annotation")
dev.off()
































































 

