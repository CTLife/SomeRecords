###########################################################################
## 1. Read the raw coverage files.
###########################################################################
myOutDir <- "100-diffMe-1-61vs64"

myFileLists <- list(
"1-Coverage-CpG/61_NC-BS2-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/61_NC-BS2-D-Girl-merge_Rep1.bismark.cov",

"1-Coverage-CpG/64_ART-BS18-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/64_ART-BS18-D-Girl_Rep1.bismark.cov"
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
               mincov = 1,       ## >= n
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
filtered.myobj_1one = filterByCoverage(myobj,  lo.count=5,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)

sink( file=paste(myOutDir_sub1, "4B-all-files-5reads.txt", sep="/") )
    print( filtered.myobj_1one )
sink()



filtered.myobj_2two = filterByCoverage(myobj,  lo.count=5,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.9)

sink( file=paste(myOutDir_sub1, "4C-all-files-5reads-rmHigh.txt", sep="/") )
print( filtered.myobj_2two )
sink()



filtered.myobj_3three = filterByCoverage(myobj,  lo.count=10,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)

sink( file=paste(myOutDir_sub1, "4D-all-files-10reads.txt", sep="/") )
print( filtered.myobj_3three )
sink()



filtered.myobj_4four = filterByCoverage(myobj,  lo.count=10,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.9)

sink( file=paste(myOutDir_sub1, "4E-all-files-10reads-rmHigh.txt", sep="/") )
print( filtered.myobj_4four )
sink()





sink( file=paste(myOutDir_sub1, "5-dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print( "######################" )
  print(   myFileLists[[i]]  )
  print(   dim(myobj[[i]])  )
  print(   dim(filtered.myobj_1one[[i]])  )
  print(   dim(filtered.myobj_2two[[i]])  )
  print(   dim(filtered.myobj_3three[[i]])  )
  print(   dim(filtered.myobj_4four[[i]])  )
}
sink()



## Merging samples
## Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. 
## This provides better coverage, but only advised when looking at CpG  
## methylation (for CpH methylation this will cause wrong results in subsequent analyses).

sink(file=paste(myOutDir_sub1, "6-log-merged-overlapSites.txt", sep="/") )

meth_1one = unite( filtered.myobj_1one , destrand=FALSE  )
head( meth_1one )
dim( meth_1one )

meth_2two = unite(filtered.myobj_2two, destrand=FALSE  )
head(meth_2two)
dim(meth_2two)

meth_3three = unite(filtered.myobj_3three, destrand=FALSE  )
head(meth_3three)
dim(meth_3three)

meth_4four = unite(filtered.myobj_4four, destrand=FALSE  )
head(meth_4four)
dim(meth_4four)

sink()




sink(file=paste(myOutDir_sub1, "7-dimensions-merged-overlap.txt", sep="/") )
    print( dim(meth_1one) )
    print( dim(meth_2two) )
    print( dim(meth_3three) )
    print( dim(meth_4four) )
sink()





sink(file=paste(myOutDir_sub1, "8-log-percMethylation.txt", sep="/") )

mat_1one = percMethylation(meth_1one)
head(mat_1one)
dim(mat_1one)

mat_2two = percMethylation(meth_2two)
head(mat_2two)
dim(mat_2two)

mat_3three = percMethylation(meth_3three)
head(mat_3three)
dim(mat_3three)

mat_4four = percMethylation(meth_4four)
head(mat_4four)
dim(mat_4four)

sink()



sink(file=paste(myOutDir_sub1, "9-dimensions-methylationMatrix.txt", sep="/") )
  print( dim(mat_1one) )
  print( dim(mat_2two) )
  print( dim(mat_3three) )
  print( dim(mat_4four) )
sink()




write.table(mat_1one , 
            file = paste(myOutDir_sub1,"10B_mathylationLevel_1one_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(mat_2two , 
            file = paste(myOutDir_sub1,"10C_mathylationLevel_2two_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three , 
            file = paste(myOutDir_sub1,"10D_mathylationLevel_3three_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four , 
            file = paste(myOutDir_sub1,"10E_mathylationLevel_4four_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

###########################################################################
###########################################################################











######################################################################################################################################################
######################################################################################################################################################
myOutDir_1one = paste(myOutDir, "/2-Cov-5reads-noRemoveHigh",  sep="");
if( ! file.exists(myOutDir_1one) ) { dir.create(myOutDir_1one, recursive = TRUE) }

pdf( file=paste(myOutDir_1one, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_1one[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_1one[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_1one[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_1one[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_1one, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_1one, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one, screeplot=TRUE)
PCASamples(meth_1one)
dev.off()




##################
myDiff_1one = calculateDiffMeth(meth_1one, num.cores=8)

myDiff25p.hypo_1one  = getMethylDiff(myDiff_1one, difference=20, qvalue=0.05, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one = getMethylDiff(myDiff_1one, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_1one       = getMethylDiff(myDiff_1one, difference=20, qvalue=0.05)


write.table(myDiff_1one , file = paste(myOutDir_1one,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one , file = paste(myOutDir_1one,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one , file = paste(myOutDir_1one,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one , file = paste(myOutDir_1one,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_1one, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_1one, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()


##########################
## Annotating differentially methylated bases or regions
myRefSeqGenes = "/home/yp/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands  = "/home/yp/AnnotationBED/hg38_CpG_islands.bed"
myRepeats     = "/home/yp/AnnotationBED/hg38_Repeats_rmsk.bed"
myImprintedRegions1 = "/home/yp/AnnotationBED/67.Regions.PlosGenetics.ImprintedGenes.hg38.bed"
myImprintedRegions2 = "/home/yp/AnnotationBED/75Regions.GR.hg38.bed"
myImprintedRegions3 = "/home/yp/AnnotationBED/369Regions.GR.hg38.bed"
myImprintedRegions4 = "/home/yp/AnnotationBED/merge1.imprintedRegions.hg38.bed"
myImprintedRegions5 = "/home/yp/AnnotationBED/merge2.imprintedRegions.hg38.bed"

gene.obj     = readTranscriptFeatures(myRefSeqGenes)
cpg.obj      = readFeatureFlank(myCpGIslands, flank=10000, feature.flank.name=c("CpGi", "shores"))
myrepeat.obj   = readFeatureFlank(myRepeats,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj = readFeatureFlank(myImprintedRegions1, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj = readFeatureFlank(myImprintedRegions2, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj = readFeatureFlank(myImprintedRegions3, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj = readFeatureFlank(myImprintedRegions4, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj = readFeatureFlank(myImprintedRegions5, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################



##
diffGeneAnn_1one_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1one,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_hypo,  percentage=TRUE)
print(diffGeneAnn_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one,"GRanges"),
                                    cpg.obj$CpGi,  cpg.obj$shores,
                                    feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_hypo,  percentage=TRUE)
print(diffCpGann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one,"GRanges"),
                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_hypo,  percentage=TRUE)
print(diffrepeatann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one, "GRanges"),
                                                       feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_hypo,  percentage=TRUE)
print(diffImprinted1Ann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one, "GRanges"),
                                                       feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_hypo,  percentage=TRUE)
print(diffImprinted2Ann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one, "GRanges"),
                                                       feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_hypo,  percentage=TRUE)
print(diffImprinted3Ann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one, "GRanges"),
                                                       feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_hypo,  percentage=TRUE)
print(diffImprinted4Ann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one, "GRanges"),
                                                       feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_hypo,  percentage=TRUE)
print(diffImprinted5Ann_1one_hypo)
sink()

pdf( file=paste(myOutDir_1one, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()






### Tiling windows analysis
tiles_1one=tileMethylCounts(filtered.myobj_1one, win.size=1000, step.size=1000)
head(tiles_1one[[1]], 3)

meth_1one_region = unite(tiles_1one, destrand=FALSE  )
head(meth_1one_region)
dim(meth_1one_region)

mat_1one_region = percMethylation(meth_1one_region)
head(mat_1one_region)
dim(mat_1one_region)

write.table(meth_1one_region , 
            file = paste(myOutDir_1one,"10A_meth_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_region , 
            file = paste(myOutDir_1one,"10B_mat_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_1one_region = calculateDiffMeth(meth_1one_region, num.cores=8)

myDiff25p.hypo_1one_region  = getMethylDiff(myDiff_1one_region, difference=20, qvalue=0.05, type="hypo")   ## less enrich in ART
myDiff25p.hyper_1one_region = getMethylDiff(myDiff_1one_region, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_1one_region       = getMethylDiff(myDiff_1one_region, difference=20, qvalue=0.05)


write.table(myDiff_1one_region , file = paste(myOutDir_1one,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one_region , file = paste(myOutDir_1one,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one_region , file = paste(myOutDir_1one,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one_region , file = paste(myOutDir_1one,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_1one, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one_region,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_1one, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one_region,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()


## Annotating differentially methylated bases or regions

##
diffGeneAnn_1one_hypo_region = annotateWithGeneParts(as(myDiff25p.hypo_1one_region, "GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_hypo_region,  percentage=TRUE)
print(diffGeneAnn_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_region,"GRanges"),
                                                cpg.obj$CpGi,  cpg.obj$shores,
                                                feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_hypo_region,  percentage=TRUE)
print(diffCpGann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_region,"GRanges"),
                                                   myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                   feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_hypo_region,  percentage=TRUE)
print(diffrepeatann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_region, "GRanges"),
                                                       feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "14A-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_hypo_region,  percentage=TRUE)
print(diffImprinted1Ann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "14B-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_region, "GRanges"),
                                                       feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "14C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_hypo_region,  percentage=TRUE)
print(diffImprinted2Ann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "14D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_region, "GRanges"),
                                                       feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "14E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_hypo_region,  percentage=TRUE)
print(diffImprinted3Ann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "14F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_region, "GRanges"),
                                                       feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "14G-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_hypo_region,  percentage=TRUE)
print(diffImprinted4Ann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "14H-distribution-onCpGs-hypo_region_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_region, "GRanges"),
                                                       feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "14I-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_hypo_region,  percentage=TRUE)
print(diffImprinted5Ann_1one_hypo_region)
sink()

pdf( file=paste(myOutDir_1one, "14J-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()
############################################################################################################################













######################################################################################################################################################
######################################################################################################################################################
myOutDir_2two = paste(myOutDir, "/3-Cov-5reads-RemoveHigh",  sep="");
if( ! file.exists(myOutDir_2two) ) { dir.create(myOutDir_2two, recursive = TRUE) }

pdf( file=paste(myOutDir_2two, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_2two[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_2two[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_2two[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_2two[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_2two, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two, screeplot=TRUE)
PCASamples(meth_2two)
dev.off()




##################
myDiff_2two = calculateDiffMeth(meth_2two, num.cores=8)

myDiff25p.hypo_2two  = getMethylDiff(myDiff_2two, difference=20, qvalue=0.05, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two = getMethylDiff(myDiff_2two, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_2two       = getMethylDiff(myDiff_2two, difference=20, qvalue=0.05)


write.table(myDiff_2two , file = paste(myOutDir_2two,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two , file = paste(myOutDir_2two,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two , file = paste(myOutDir_2two,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two , file = paste(myOutDir_2two,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_2two, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_2two, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()




##
diffGeneAnn_2two_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two,"GRanges"),
                                                cpg.obj$CpGi,  cpg.obj$shores,
                                                feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_hypo,  percentage=TRUE)
print(diffCpGann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two,"GRanges"),
                                                   myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                   feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_hypo,  percentage=TRUE)
print(diffrepeatann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two, "GRanges"),
                                                       feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two, "GRanges"),
                                                       feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two, "GRanges"),
                                                       feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two, "GRanges"),
                                                       feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two, "GRanges"),
                                                       feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                       feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_hypo)
sink()

pdf( file=paste(myOutDir_2two, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()






### Tiling windows analysis
tiles_2two=tileMethylCounts(filtered.myobj_2two, win.size=1000, step.size=1000)
head(tiles_2two[[1]], 3)

meth_2two_region = unite(tiles_2two, destrand=FALSE  )
head(meth_2two_region)
dim(meth_2two_region)

mat_2two_region = percMethylation(meth_2two_region)
head(mat_2two_region)
dim(mat_2two_region)

write.table(meth_2two_region , 
            file = paste(myOutDir_2two,"10A_meth_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_region , 
            file = paste(myOutDir_2two,"10B_mat_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_2two_region = calculateDiffMeth(meth_2two_region, num.cores=8)

myDiff25p.hypo_2two_region  = getMethylDiff(myDiff_2two_region, difference=20, qvalue=0.05, type="hypo")   ## less enrich in ART
myDiff25p.hyper_2two_region = getMethylDiff(myDiff_2two_region, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_2two_region       = getMethylDiff(myDiff_2two_region, difference=20, qvalue=0.05)


write.table(myDiff_2two_region , file = paste(myOutDir_2two,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_region , file = paste(myOutDir_2two,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_region , file = paste(myOutDir_2two,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_region , file = paste(myOutDir_2two,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_2two, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_region,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_2two, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_region,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()


## Annotating differentially methylated bases or regions

##
diffGeneAnn_2two_hypo_region = annotateWithGeneParts(as(myDiff25p.hypo_2two_region, "GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_hypo_region,  percentage=TRUE)
print(diffGeneAnn_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_region,"GRanges"),
                                                       cpg.obj$CpGi,  cpg.obj$shores,
                                                       feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_hypo_region,  percentage=TRUE)
print(diffCpGann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_region,"GRanges"),
                                                          myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                          feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_hypo_region,  percentage=TRUE)
print(diffrepeatann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_region, "GRanges"),
                                                              feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "14A-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_hypo_region,  percentage=TRUE)
print(diffImprinted1Ann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "14B-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_region, "GRanges"),
                                                              feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "14C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_hypo_region,  percentage=TRUE)
print(diffImprinted2Ann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "14D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_region, "GRanges"),
                                                              feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "14E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_hypo_region,  percentage=TRUE)
print(diffImprinted3Ann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "14F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_region, "GRanges"),
                                                              feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "14G-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_hypo_region,  percentage=TRUE)
print(diffImprinted4Ann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "14H-distribution-onCpGs-hypo_region_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_region, "GRanges"),
                                                              feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "14I-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_hypo_region,  percentage=TRUE)
print(diffImprinted5Ann_2two_hypo_region)
sink()

pdf( file=paste(myOutDir_2two, "14J-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()
############################################################################################################################
















######################################################################################################################################################
######################################################################################################################################################
myOutDir_3three = paste(myOutDir, "/4-Cov-10reads-noRemoveHigh",  sep="");
if( ! file.exists(myOutDir_3three) ) { dir.create(myOutDir_3three, recursive = TRUE) }

pdf( file=paste(myOutDir_3three, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_3three[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_3three[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_3three[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_3three[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_3three, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three, screeplot=TRUE)
PCASamples(meth_3three)
dev.off()




##################
myDiff_3three = calculateDiffMeth(meth_3three, num.cores=8)

myDiff25p.hypo_3three  = getMethylDiff(myDiff_3three, difference=20, qvalue=0.05, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three = getMethylDiff(myDiff_3three, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_3three       = getMethylDiff(myDiff_3three, difference=20, qvalue=0.05)


write.table(myDiff_3three , file = paste(myOutDir_3three,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three , file = paste(myOutDir_3three,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three , file = paste(myOutDir_3three,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three , file = paste(myOutDir_3three,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_3three, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_3three, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()




##
diffGeneAnn_3three_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three,"GRanges"),
                                                  cpg.obj$CpGi,  cpg.obj$shores,
                                                  feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_hypo,  percentage=TRUE)
print(diffCpGann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three,"GRanges"),
                                                     myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                     feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_hypo,  percentage=TRUE)
print(diffrepeatann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three, "GRanges"),
                                                         feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three, "GRanges"),
                                                         feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three, "GRanges"),
                                                         feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three, "GRanges"),
                                                         feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three, "GRanges"),
                                                         feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_hypo)
sink()

pdf( file=paste(myOutDir_3three, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()






### Tiling windows analysis
tiles_3three=tileMethylCounts(filtered.myobj_3three, win.size=1000, step.size=1000)
head(tiles_3three[[1]], 3)

meth_3three_region = unite(tiles_3three, destrand=FALSE  )
head(meth_3three_region)
dim(meth_3three_region)

mat_3three_region = percMethylation(meth_3three_region)
head(mat_3three_region)
dim(mat_3three_region)

write.table(meth_3three_region , 
            file = paste(myOutDir_3three,"10A_meth_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_region , 
            file = paste(myOutDir_3three,"10B_mat_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_3three_region = calculateDiffMeth(meth_3three_region, num.cores=8)

myDiff25p.hypo_3three_region  = getMethylDiff(myDiff_3three_region, difference=20, qvalue=0.05, type="hypo")   ## less enrich in ART
myDiff25p.hyper_3three_region = getMethylDiff(myDiff_3three_region, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_3three_region       = getMethylDiff(myDiff_3three_region, difference=20, qvalue=0.05)


write.table(myDiff_3three_region , file = paste(myOutDir_3three,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_region , file = paste(myOutDir_3three,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_region , file = paste(myOutDir_3three,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_region , file = paste(myOutDir_3three,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_3three, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_region,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_3three, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_region,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()


## Annotating differentially methylated bases or regions

##
diffGeneAnn_3three_hypo_region = annotateWithGeneParts(as(myDiff25p.hypo_3three_region, "GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_hypo_region,  percentage=TRUE)
print(diffGeneAnn_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_region,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_hypo_region,  percentage=TRUE)
print(diffCpGann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_region,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_hypo_region,  percentage=TRUE)
print(diffrepeatann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_region, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "14A-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_hypo_region,  percentage=TRUE)
print(diffImprinted1Ann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "14B-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_region, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "14C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_hypo_region,  percentage=TRUE)
print(diffImprinted2Ann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "14D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_region, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "14E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_hypo_region,  percentage=TRUE)
print(diffImprinted3Ann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "14F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_region, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "14G-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_hypo_region,  percentage=TRUE)
print(diffImprinted4Ann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "14H-distribution-onCpGs-hypo_region_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_region, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "14I-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_hypo_region,  percentage=TRUE)
print(diffImprinted5Ann_3three_hypo_region)
sink()

pdf( file=paste(myOutDir_3three, "14J-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()
############################################################################################################################






















######################################################################################################################################################
######################################################################################################################################################
myOutDir_4four = paste(myOutDir, "/5-Cov-10reads-RemoveHigh",  sep="");
if( ! file.exists(myOutDir_4four) ) { dir.create(myOutDir_4four, recursive = TRUE) }

pdf( file=paste(myOutDir_4four, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(filtered.myobj_4four[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( filtered.myobj_4four[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(filtered.myobj_4four[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( filtered.myobj_4four[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_4four, "3-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four, "4-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four, screeplot=TRUE)
PCASamples(meth_4four)
dev.off()




##################
myDiff_4four = calculateDiffMeth(meth_4four, num.cores=8)

myDiff25p.hypo_4four  = getMethylDiff(myDiff_4four, difference=20, qvalue=0.05, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four = getMethylDiff(myDiff_4four, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_4four       = getMethylDiff(myDiff_4four, difference=20, qvalue=0.05)


write.table(myDiff_4four , file = paste(myOutDir_4four,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four , file = paste(myOutDir_4four,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four , file = paste(myOutDir_4four,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four , file = paste(myOutDir_4four,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_4four, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_4four, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()




##
diffGeneAnn_4four_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four,"GRanges"),
                                                 cpg.obj$CpGi,  cpg.obj$shores,
                                                 feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_hypo,  percentage=TRUE)
print(diffCpGann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four,"GRanges"),
                                                    myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                    feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_hypo,  percentage=TRUE)
print(diffrepeatann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four, "GRanges"),
                                                        feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four, "GRanges"),
                                                        feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four, "GRanges"),
                                                        feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four, "GRanges"),
                                                        feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four, "GRanges"),
                                                        feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_hypo)
sink()

pdf( file=paste(myOutDir_4four, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()






### Tiling windows analysis
tiles_4four=tileMethylCounts(filtered.myobj_4four, win.size=1000, step.size=1000)
head(tiles_4four[[1]], 3)

meth_4four_region = unite(tiles_4four, destrand=FALSE  )
head(meth_4four_region)
dim(meth_4four_region)

mat_4four_region = percMethylation(meth_4four_region)
head(mat_4four_region)
dim(mat_4four_region)

write.table(meth_4four_region , 
            file = paste(myOutDir_4four,"10A_meth_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_region , 
            file = paste(myOutDir_4four,"10B_mat_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



myDiff_4four_region = calculateDiffMeth(meth_4four_region, num.cores=8)

myDiff25p.hypo_4four_region  = getMethylDiff(myDiff_4four_region, difference=20, qvalue=0.05, type="hypo")   ## less enrich in ART
myDiff25p.hyper_4four_region = getMethylDiff(myDiff_4four_region, difference=20, qvalue=0.05, type="hyper")  ## more enrich in ART
myDiff25p_4four_region       = getMethylDiff(myDiff_4four_region, difference=20, qvalue=0.05)


write.table(myDiff_4four_region , file = paste(myOutDir_4four,"11A_diffMe-allsites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_region , file = paste(myOutDir_4four,"11B_diffMe-hypo_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_region , file = paste(myOutDir_4four,"11C_diffMe-hyper_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_region , file = paste(myOutDir_4four,"11D_AlldiffMesites_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



sink( file=paste(myOutDir_4four, "12A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_region,  plot=FALSE,  qvalue.cutoff=0.05, meth.cutoff=20)
sink()


pdf( file=paste(myOutDir_4four, "12B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_region,  plot=TRUE,  qvalue.cutoff=0.05, meth.cutoff=20)
dev.off()


## Annotating differentially methylated bases or regions

##
diffGeneAnn_4four_hypo_region = annotateWithGeneParts(as(myDiff25p.hypo_4four_region, "GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four, "13A-distribution-onGenes-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_hypo_region,  percentage=TRUE)
print(diffGeneAnn_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "13B-distribution-onGenes-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_hypo_region,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_region,"GRanges"),
                                                        cpg.obj$CpGi,  cpg.obj$shores,
                                                        feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four, "13C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_hypo_region,  percentage=TRUE)
print(diffCpGann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "13D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_hypo_region = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_region,"GRanges"),
                                                           myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                           feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four, "13E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_hypo_region,  percentage=TRUE)
print(diffrepeatann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "13F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_region, "GRanges"),
                                                               feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "14A-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_hypo_region,  percentage=TRUE)
print(diffImprinted1Ann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "14B-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_region, "GRanges"),
                                                               feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "14C-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_hypo_region,  percentage=TRUE)
print(diffImprinted2Ann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "14D-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_region, "GRanges"),
                                                               feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "14E-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_hypo_region,  percentage=TRUE)
print(diffImprinted3Ann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "14F-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_region, "GRanges"),
                                                               feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "14G-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_hypo_region,  percentage=TRUE)
print(diffImprinted4Ann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "14H-distribution-onCpGs-hypo_region_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_hypo_region = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_region, "GRanges"),
                                                               feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "14I-distribution-on-hypo_region.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_hypo_region,  percentage=TRUE)
print(diffImprinted5Ann_4four_hypo_region)
sink()

pdf( file=paste(myOutDir_4four, "14J-distribution-onCpGs-hypo_region.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_hypo_region, precedence=TRUE, main="differential methylation annotation")
dev.off()
############################################################################################################################















































