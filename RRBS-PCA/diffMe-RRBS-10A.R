###########################################################################
## 1. Read the raw coverage files.
###########################################################################
myOutDir <- "10A_diffMe-FemaleTwin_NC-vs-IVFfresh"

myFileLists <- list(
"1-Coverage-CpG/61_NC-BS2-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/61_NC-BS2-D-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/62_NC-BS20-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/62_NC-BS20-D-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/63_NC-E8-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/63_NC-E8-D-Girl-merge_Rep1.bismark.cov",

"1-Coverage-CpG/64_ART-BS18-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/64_ART-BS18-D-Girl_Rep1.bismark.cov",
"1-Coverage-CpG/65_ART-BS29-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/65_ART-BS29-D-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/66_ART-E29-C-Girl-merge_Rep1.bismark.cov",
"1-Coverage-CpG/66_ART-E29-D-Girl-merge_Rep1.bismark.cov"
)

mySampleID <- list( "NC1", "NC2",  "NC3", "NC4", "NC5", "NC6", "IVF-fresh1", "IVF-fresh2" 
                    , "IVF-fresh3", "IVF-fresh4" , "IVF-fresh5", "IVF-fresh6"  )
myTreatment <- c( 0, 0,  0, 0,  0, 0, 1, 1, 1, 1, 1, 1)       ## This option will determine the result of unite.




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

meth_1one = unite( filtered.myobj_1one , destrand=FALSE,  min.per.group = 5L   )
head( meth_1one )
dim( meth_1one )

meth_2two = unite(filtered.myobj_2two, destrand=FALSE,  min.per.group = 5L   )
head(meth_2two)
dim(meth_2two)

meth_3three = unite(filtered.myobj_3three, destrand=FALSE,  min.per.group = 5L   )
head(meth_3three)
dim(meth_3three)

meth_4four = unite(filtered.myobj_4four, destrand=FALSE,  min.per.group = 5L   )
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
sink()


pdf( file=paste(myOutDir_1one, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_1one, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one, screeplot=TRUE)
PCASamples(meth_1one)
dev.off()



##################
myDiff_1one = calculateDiffMeth(meth_1one, num.cores=16)
dim(myDiff_1one)
names(myDiff_1one)

myQvalue_1one = myDiff_1one$qvalue
myMethDi_1one = myDiff_1one$meth.diff

pdf(paste(myOutDir_1one,"4A_qvalue_distribution.pdf",  sep="/"))
    hist(myQvalue_1one, nclass=100, xlim=c(0, 1), freq=FALSE)
    hist(myQvalue_1one[myQvalue_1one<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_1one,  "4B_methDiff_distribution.pdf",  sep="/"))
    hist(myMethDi_1one, nclass=100, xlim=c(0, 100), freq=FALSE)
    hist(myMethDi_1one[myMethDi_1one>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
    hist(myMethDi_1one[myMethDi_1one>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
    hist(myMethDi_1one[myMethDi_1one<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_1one <- -log10(myQvalue_1one)
myColor_1one[myColor_1one >  2]  = "yes"
myColor_1one[myColor_1one <= 2]  = "no"
DataFrame_Local5_1one <- data.frame(myx1 = -log10(myQvalue_1one), myy1 = myMethDi_1one,   mycolor1 = myColor_1one )
  
ggplot(DataFrame_Local5_1one, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
                geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
                xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_1one,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_1one,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_1one,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()

 



myDiff25p.hypo_1one  = getMethylDiff(myDiff_1one, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one = getMethylDiff(myDiff_1one, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_1one       = getMethylDiff(myDiff_1one, difference=10, qvalue=0.01)
myDiffTemp_1one      = getMethylDiff(myDiff_1one, difference=0,  qvalue=0.05)

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

write.table(myDiffTemp_1one , file = paste(myOutDir_1one,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_1one, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_1one, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
myrepeat.obj = readFeatureFlank(myRepeats,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj = readFeatureFlank(myImprintedRegions1, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj = readFeatureFlank(myImprintedRegions2, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj = readFeatureFlank(myImprintedRegions3, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj = readFeatureFlank(myImprintedRegions4, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj = readFeatureFlank(myImprintedRegions5, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################



########################## annotation for hypo sites.
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
####################






########################## annotation for hyper sites.
diffGeneAnn_1one_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1one,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_hyper,  percentage=TRUE)
print(diffGeneAnn_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one,"GRanges"),
                                                 cpg.obj$CpGi,  cpg.obj$shores,
                                                 feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_hyper,  percentage=TRUE)
print(diffCpGann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one,"GRanges"),
                                                    myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                    feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_hyper,  percentage=TRUE)
print(diffrepeatann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one, "GRanges"),
                                                        feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_hyper,  percentage=TRUE)
print(diffImprinted1Ann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one, "GRanges"),
                                                        feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_hyper,  percentage=TRUE)
print(diffImprinted2Ann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one, "GRanges"),
                                                        feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_hyper,  percentage=TRUE)
print(diffImprinted3Ann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one, "GRanges"),
                                                        feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_hyper,  percentage=TRUE)
print(diffImprinted4Ann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one, "GRanges"),
                                                        feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_hyper,  percentage=TRUE)
print(diffImprinted5Ann_1one_hyper)
sink()

pdf( file=paste(myOutDir_1one, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_1one_sub1_10kb = paste(myOutDir_1one, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_1one_sub1_10kb) ) { dir.create(myOutDir_1one_sub1_10kb, recursive = TRUE) }


tiles_1one_sub1_10kb = tileMethylCounts(filtered.myobj_1one, win.size=10000, step.size=10000)
head(tiles_1one_sub1_10kb[[1]], 3)

meth_1one_sub1_10kb = unite(tiles_1one_sub1_10kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_1one_sub1_10kb)
dim(meth_1one_sub1_10kb)

mat_1one_sub1_10kb = percMethylation(meth_1one_sub1_10kb)
head(mat_1one_sub1_10kb)
dim(mat_1one_sub1_10kb)

write.table(meth_1one_sub1_10kb , 
            file = paste(myOutDir_1one_sub1_10kb,   "0A_meth_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_sub1_10kb , 
            file = paste(myOutDir_1one_sub1_10kb,   "0B_mat_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_1one_sub1_10kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_1one_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one_sub1_10kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_1one_sub1_10kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one_sub1_10kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_1one_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one_sub1_10kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_1one_sub1_10kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_1one_sub1_10kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_1one_sub1_10kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub1_10kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub1_10kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_1one_sub1_10kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_1one_sub1_10kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_1one_sub1_10kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub1_10kb, screeplot=TRUE)
PCASamples(meth_1one_sub1_10kb)
dev.off()



##################
myDiff_1one_sub1_10kb = calculateDiffMeth(meth_1one_sub1_10kb, num.cores=16)
dim(myDiff_1one_sub1_10kb)
names(myDiff_1one_sub1_10kb)

myQvalue_1one_sub1_10kb = myDiff_1one_sub1_10kb$qvalue
myMethDi_1one_sub1_10kb = myDiff_1one_sub1_10kb$meth.diff

pdf(paste(myOutDir_1one_sub1_10kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_1one_sub1_10kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_1one_sub1_10kb[myQvalue_1one_sub1_10kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_1one_sub1_10kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_1one_sub1_10kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_1one_sub1_10kb[myMethDi_1one_sub1_10kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_1one_sub1_10kb[myMethDi_1one_sub1_10kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_1one_sub1_10kb[myMethDi_1one_sub1_10kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_1one_sub1_10kb <- -log10(myQvalue_1one_sub1_10kb)
myColor_1one_sub1_10kb[myColor_1one_sub1_10kb >  2]  = "yes"
myColor_1one_sub1_10kb[myColor_1one_sub1_10kb <= 2]  = "no"
DataFrame_Local5_1one_sub1_10kb <- data.frame(myx1 = -log10(myQvalue_1one_sub1_10kb), myy1 = myMethDi_1one_sub1_10kb,   mycolor1 = myColor_1one_sub1_10kb )

ggplot(DataFrame_Local5_1one_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_1one_sub1_10kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_1one_sub1_10kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_1one_sub1_10kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_1one_sub1_10kb  = getMethylDiff(myDiff_1one_sub1_10kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one_sub1_10kb = getMethylDiff(myDiff_1one_sub1_10kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_1one_sub1_10kb       = getMethylDiff(myDiff_1one_sub1_10kb, difference=10, qvalue=0.01)
myDiffTemp_1one_sub1_10kb      = getMethylDiff(myDiff_1one_sub1_10kb, difference=0,  qvalue=0.05)

write.table(myDiff_1one_sub1_10kb , file = paste(myOutDir_1one_sub1_10kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one_sub1_10kb , file = paste(myOutDir_1one_sub1_10kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one_sub1_10kb , file = paste(myOutDir_1one_sub1_10kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one_sub1_10kb , file = paste(myOutDir_1one_sub1_10kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_1one_sub1_10kb , file = paste(myOutDir_1one_sub1_10kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_1one_sub1_10kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one_sub1_10kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_1one_sub1_10kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one_sub1_10kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()


 

########################## annotation for hypo sites.
diffGeneAnn_1one_sub1_10kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1one_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub1_10kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffGeneAnn_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub1_10kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub1_10kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffCpGann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub1_10kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffrepeatann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub1_10kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub1_10kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub1_10kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub1_10kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub1_10kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_1one_sub1_10kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1one_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub1_10kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffGeneAnn_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub1_10kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub1_10kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffCpGann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub1_10kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffrepeatann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub1_10kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub1_10kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub1_10kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub1_10kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub1_10kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub1_10kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub1_10kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()







### Tiling windows analysis
##################
myOutDir_1one_sub2_5kb = paste(myOutDir_1one, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_1one_sub2_5kb) ) { dir.create(myOutDir_1one_sub2_5kb, recursive = TRUE) }


tiles_1one_sub2_5kb = tileMethylCounts(filtered.myobj_1one, win.size=5000, step.size=5000)
head(tiles_1one_sub2_5kb[[1]], 3)

meth_1one_sub2_5kb = unite(tiles_1one_sub2_5kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_1one_sub2_5kb)
dim(meth_1one_sub2_5kb)

mat_1one_sub2_5kb = percMethylation(meth_1one_sub2_5kb)
head(mat_1one_sub2_5kb)
dim(mat_1one_sub2_5kb)

write.table(meth_1one_sub2_5kb , 
            file = paste(myOutDir_1one_sub2_5kb,   "0A_meth_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_sub2_5kb , 
            file = paste(myOutDir_1one_sub2_5kb,   "0B_mat_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_1one_sub2_5kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_1one_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one_sub2_5kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_1one_sub2_5kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one_sub2_5kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_1one_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one_sub2_5kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_1one_sub2_5kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_1one_sub2_5kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_1one_sub2_5kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub2_5kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub2_5kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_1one_sub2_5kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_1one_sub2_5kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_1one_sub2_5kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub2_5kb, screeplot=TRUE)
PCASamples(meth_1one_sub2_5kb)
dev.off()



##################
myDiff_1one_sub2_5kb = calculateDiffMeth(meth_1one_sub2_5kb, num.cores=16)
dim(myDiff_1one_sub2_5kb)
names(myDiff_1one_sub2_5kb)

myQvalue_1one_sub2_5kb = myDiff_1one_sub2_5kb$qvalue
myMethDi_1one_sub2_5kb = myDiff_1one_sub2_5kb$meth.diff

pdf(paste(myOutDir_1one_sub2_5kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_1one_sub2_5kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_1one_sub2_5kb[myQvalue_1one_sub2_5kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_1one_sub2_5kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_1one_sub2_5kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_1one_sub2_5kb[myMethDi_1one_sub2_5kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_1one_sub2_5kb[myMethDi_1one_sub2_5kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_1one_sub2_5kb[myMethDi_1one_sub2_5kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_1one_sub2_5kb <- -log10(myQvalue_1one_sub2_5kb)
myColor_1one_sub2_5kb[myColor_1one_sub2_5kb >  2]  = "yes"
myColor_1one_sub2_5kb[myColor_1one_sub2_5kb <= 2]  = "no"
DataFrame_Local5_1one_sub2_5kb <- data.frame(myx1 = -log10(myQvalue_1one_sub2_5kb), myy1 = myMethDi_1one_sub2_5kb,   mycolor1 = myColor_1one_sub2_5kb )

ggplot(DataFrame_Local5_1one_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_1one_sub2_5kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_1one_sub2_5kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_1one_sub2_5kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_1one_sub2_5kb  = getMethylDiff(myDiff_1one_sub2_5kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one_sub2_5kb = getMethylDiff(myDiff_1one_sub2_5kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_1one_sub2_5kb       = getMethylDiff(myDiff_1one_sub2_5kb, difference=10, qvalue=0.01)
myDiffTemp_1one_sub2_5kb      = getMethylDiff(myDiff_1one_sub2_5kb, difference=0,  qvalue=0.05)

write.table(myDiff_1one_sub2_5kb , file = paste(myOutDir_1one_sub2_5kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one_sub2_5kb , file = paste(myOutDir_1one_sub2_5kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one_sub2_5kb , file = paste(myOutDir_1one_sub2_5kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one_sub2_5kb , file = paste(myOutDir_1one_sub2_5kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_1one_sub2_5kb , file = paste(myOutDir_1one_sub2_5kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_1one_sub2_5kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one_sub2_5kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_1one_sub2_5kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one_sub2_5kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_1one_sub2_5kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1one_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub2_5kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffGeneAnn_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub2_5kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub2_5kb,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffCpGann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub2_5kb,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffrepeatann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub2_5kb, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub2_5kb, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub2_5kb, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub2_5kb, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub2_5kb, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_1one_sub2_5kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1one_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub2_5kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffGeneAnn_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub2_5kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub2_5kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffCpGann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub2_5kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffrepeatann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub2_5kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub2_5kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub2_5kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub2_5kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub2_5kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub2_5kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub2_5kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()










### Tiling windows analysis
##################
myOutDir_1one_sub3_1kb = paste(myOutDir_1one, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_1one_sub3_1kb) ) { dir.create(myOutDir_1one_sub3_1kb, recursive = TRUE) }


tiles_1one_sub3_1kb = tileMethylCounts(filtered.myobj_1one, win.size=1000, step.size=1000)
head(tiles_1one_sub3_1kb[[1]], 3)

meth_1one_sub3_1kb = unite(tiles_1one_sub3_1kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_1one_sub3_1kb)
dim(meth_1one_sub3_1kb)

mat_1one_sub3_1kb = percMethylation(meth_1one_sub3_1kb)
head(mat_1one_sub3_1kb)
dim(mat_1one_sub3_1kb)

write.table(meth_1one_sub3_1kb , 
            file = paste(myOutDir_1one_sub3_1kb,   "0A_meth_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_sub3_1kb , 
            file = paste(myOutDir_1one_sub3_1kb,   "0B_mat_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_1one_sub3_1kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_1one_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one_sub3_1kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_1one_sub3_1kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one_sub3_1kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_1one_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one_sub3_1kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_1one_sub3_1kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_1one_sub3_1kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_1one_sub3_1kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub3_1kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub3_1kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_1one_sub3_1kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_1one_sub3_1kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_1one_sub3_1kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub3_1kb, screeplot=TRUE)
PCASamples(meth_1one_sub3_1kb)
dev.off()



##################
myDiff_1one_sub3_1kb = calculateDiffMeth(meth_1one_sub3_1kb, num.cores=16)
dim(myDiff_1one_sub3_1kb)
names(myDiff_1one_sub3_1kb)

myQvalue_1one_sub3_1kb = myDiff_1one_sub3_1kb$qvalue
myMethDi_1one_sub3_1kb = myDiff_1one_sub3_1kb$meth.diff

pdf(paste(myOutDir_1one_sub3_1kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_1one_sub3_1kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_1one_sub3_1kb[myQvalue_1one_sub3_1kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_1one_sub3_1kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_1one_sub3_1kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_1one_sub3_1kb[myMethDi_1one_sub3_1kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_1one_sub3_1kb[myMethDi_1one_sub3_1kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_1one_sub3_1kb[myMethDi_1one_sub3_1kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_1one_sub3_1kb <- -log10(myQvalue_1one_sub3_1kb)
myColor_1one_sub3_1kb[myColor_1one_sub3_1kb >  2]  = "yes"
myColor_1one_sub3_1kb[myColor_1one_sub3_1kb <= 2]  = "no"
DataFrame_Local5_1one_sub3_1kb <- data.frame(myx1 = -log10(myQvalue_1one_sub3_1kb), myy1 = myMethDi_1one_sub3_1kb,   mycolor1 = myColor_1one_sub3_1kb )

ggplot(DataFrame_Local5_1one_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_1one_sub3_1kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_1one_sub3_1kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_1one_sub3_1kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_1one_sub3_1kb  = getMethylDiff(myDiff_1one_sub3_1kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one_sub3_1kb = getMethylDiff(myDiff_1one_sub3_1kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_1one_sub3_1kb       = getMethylDiff(myDiff_1one_sub3_1kb, difference=10, qvalue=0.01)
myDiffTemp_1one_sub3_1kb      = getMethylDiff(myDiff_1one_sub3_1kb, difference=0,  qvalue=0.05)

write.table(myDiff_1one_sub3_1kb , file = paste(myOutDir_1one_sub3_1kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one_sub3_1kb , file = paste(myOutDir_1one_sub3_1kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one_sub3_1kb , file = paste(myOutDir_1one_sub3_1kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one_sub3_1kb , file = paste(myOutDir_1one_sub3_1kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_1one_sub3_1kb , file = paste(myOutDir_1one_sub3_1kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_1one_sub3_1kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one_sub3_1kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_1one_sub3_1kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one_sub3_1kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_1one_sub3_1kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1one_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub3_1kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffGeneAnn_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub3_1kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub3_1kb,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffCpGann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub3_1kb,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffrepeatann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub3_1kb, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub3_1kb, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub3_1kb, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub3_1kb, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub3_1kb, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_1one_sub3_1kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1one_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub3_1kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffGeneAnn_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub3_1kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub3_1kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffCpGann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub3_1kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffrepeatann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub3_1kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub3_1kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub3_1kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub3_1kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub3_1kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub3_1kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub3_1kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_1one_sub4_500bp = paste(myOutDir_1one, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_1one_sub4_500bp) ) { dir.create(myOutDir_1one_sub4_500bp, recursive = TRUE) }


tiles_1one_sub4_500bp = tileMethylCounts(filtered.myobj_1one, win.size=500, step.size=500)
head(tiles_1one_sub4_500bp[[1]], 3)

meth_1one_sub4_500bp = unite(tiles_1one_sub4_500bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_1one_sub4_500bp)
dim(meth_1one_sub4_500bp)

mat_1one_sub4_500bp = percMethylation(meth_1one_sub4_500bp)
head(mat_1one_sub4_500bp)
dim(mat_1one_sub4_500bp)

write.table(meth_1one_sub4_500bp , 
            file = paste(myOutDir_1one_sub4_500bp,   "0A_meth_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_sub4_500bp , 
            file = paste(myOutDir_1one_sub4_500bp,   "0B_mat_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_1one_sub4_500bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_1one_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one_sub4_500bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_1one_sub4_500bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one_sub4_500bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_1one_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one_sub4_500bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_1one_sub4_500bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_1one_sub4_500bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_1one_sub4_500bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub4_500bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub4_500bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_1one_sub4_500bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_1one_sub4_500bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_1one_sub4_500bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub4_500bp, screeplot=TRUE)
PCASamples(meth_1one_sub4_500bp)
dev.off()



##################
myDiff_1one_sub4_500bp = calculateDiffMeth(meth_1one_sub4_500bp, num.cores=16)
dim(myDiff_1one_sub4_500bp)
names(myDiff_1one_sub4_500bp)

myQvalue_1one_sub4_500bp = myDiff_1one_sub4_500bp$qvalue
myMethDi_1one_sub4_500bp = myDiff_1one_sub4_500bp$meth.diff

pdf(paste(myOutDir_1one_sub4_500bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_1one_sub4_500bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_1one_sub4_500bp[myQvalue_1one_sub4_500bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_1one_sub4_500bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_1one_sub4_500bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_1one_sub4_500bp[myMethDi_1one_sub4_500bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_1one_sub4_500bp[myMethDi_1one_sub4_500bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_1one_sub4_500bp[myMethDi_1one_sub4_500bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_1one_sub4_500bp <- -log10(myQvalue_1one_sub4_500bp)
myColor_1one_sub4_500bp[myColor_1one_sub4_500bp >  2]  = "yes"
myColor_1one_sub4_500bp[myColor_1one_sub4_500bp <= 2]  = "no"
DataFrame_Local5_1one_sub4_500bp <- data.frame(myx1 = -log10(myQvalue_1one_sub4_500bp), myy1 = myMethDi_1one_sub4_500bp,   mycolor1 = myColor_1one_sub4_500bp )

ggplot(DataFrame_Local5_1one_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_1one_sub4_500bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_1one_sub4_500bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_1one_sub4_500bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_1one_sub4_500bp  = getMethylDiff(myDiff_1one_sub4_500bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one_sub4_500bp = getMethylDiff(myDiff_1one_sub4_500bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_1one_sub4_500bp       = getMethylDiff(myDiff_1one_sub4_500bp, difference=10, qvalue=0.01)
myDiffTemp_1one_sub4_500bp      = getMethylDiff(myDiff_1one_sub4_500bp, difference=0,  qvalue=0.05)

write.table(myDiff_1one_sub4_500bp , file = paste(myOutDir_1one_sub4_500bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one_sub4_500bp , file = paste(myOutDir_1one_sub4_500bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one_sub4_500bp , file = paste(myOutDir_1one_sub4_500bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one_sub4_500bp , file = paste(myOutDir_1one_sub4_500bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_1one_sub4_500bp , file = paste(myOutDir_1one_sub4_500bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_1one_sub4_500bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one_sub4_500bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_1one_sub4_500bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one_sub4_500bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_1one_sub4_500bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1one_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub4_500bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffGeneAnn_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub4_500bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub4_500bp,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffCpGann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub4_500bp,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffrepeatann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub4_500bp, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub4_500bp, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub4_500bp, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub4_500bp, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub4_500bp, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_1one_sub4_500bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1one_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub4_500bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffGeneAnn_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub4_500bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub4_500bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffCpGann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub4_500bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffrepeatann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub4_500bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub4_500bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub4_500bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub4_500bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub4_500bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub4_500bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub4_500bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
















### Tiling windows analysis
##################
myOutDir_1one_sub5_100bp = paste(myOutDir_1one, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_1one_sub5_100bp) ) { dir.create(myOutDir_1one_sub5_100bp, recursive = TRUE) }


tiles_1one_sub5_100bp = tileMethylCounts(filtered.myobj_1one, win.size=100, step.size=100)
head(tiles_1one_sub5_100bp[[1]], 3)

meth_1one_sub5_100bp = unite(tiles_1one_sub5_100bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_1one_sub5_100bp)
dim(meth_1one_sub5_100bp)

mat_1one_sub5_100bp = percMethylation(meth_1one_sub5_100bp)
head(mat_1one_sub5_100bp)
dim(mat_1one_sub5_100bp)

write.table(meth_1one_sub5_100bp , 
            file = paste(myOutDir_1one_sub5_100bp,   "0A_meth_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_sub5_100bp , 
            file = paste(myOutDir_1one_sub5_100bp,   "0B_mat_1one_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_1one_sub5_100bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_1one_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one_sub5_100bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_1one_sub5_100bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one_sub5_100bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_1one_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one_sub5_100bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_1one_sub5_100bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_1one_sub5_100bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_1one_sub5_100bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub5_100bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_1one_sub5_100bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_1one_sub5_100bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_1one_sub5_100bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_1one_sub5_100bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_1one_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub5_100bp, screeplot=TRUE)
PCASamples(meth_1one_sub5_100bp)
dev.off()



##################
myDiff_1one_sub5_100bp = calculateDiffMeth(meth_1one_sub5_100bp, num.cores=16)
dim(myDiff_1one_sub5_100bp)
names(myDiff_1one_sub5_100bp)

myQvalue_1one_sub5_100bp = myDiff_1one_sub5_100bp$qvalue
myMethDi_1one_sub5_100bp = myDiff_1one_sub5_100bp$meth.diff

pdf(paste(myOutDir_1one_sub5_100bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_1one_sub5_100bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_1one_sub5_100bp[myQvalue_1one_sub5_100bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_1one_sub5_100bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_1one_sub5_100bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_1one_sub5_100bp[myMethDi_1one_sub5_100bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_1one_sub5_100bp[myMethDi_1one_sub5_100bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_1one_sub5_100bp[myMethDi_1one_sub5_100bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_1one_sub5_100bp <- -log10(myQvalue_1one_sub5_100bp)
myColor_1one_sub5_100bp[myColor_1one_sub5_100bp >  2]  = "yes"
myColor_1one_sub5_100bp[myColor_1one_sub5_100bp <= 2]  = "no"
DataFrame_Local5_1one_sub5_100bp <- data.frame(myx1 = -log10(myQvalue_1one_sub5_100bp), myy1 = myMethDi_1one_sub5_100bp,   mycolor1 = myColor_1one_sub5_100bp )

ggplot(DataFrame_Local5_1one_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_1one_sub5_100bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_1one_sub5_100bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_1one_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_1one_sub5_100bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_1one_sub5_100bp  = getMethylDiff(myDiff_1one_sub5_100bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_1one_sub5_100bp = getMethylDiff(myDiff_1one_sub5_100bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_1one_sub5_100bp       = getMethylDiff(myDiff_1one_sub5_100bp, difference=10, qvalue=0.01)
myDiffTemp_1one_sub5_100bp      = getMethylDiff(myDiff_1one_sub5_100bp, difference=0,  qvalue=0.05)

write.table(myDiff_1one_sub5_100bp , file = paste(myOutDir_1one_sub5_100bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_1one_sub5_100bp , file = paste(myOutDir_1one_sub5_100bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_1one_sub5_100bp , file = paste(myOutDir_1one_sub5_100bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_1one_sub5_100bp , file = paste(myOutDir_1one_sub5_100bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_1one_sub5_100bp , file = paste(myOutDir_1one_sub5_100bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_1one_sub5_100bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_1one_sub5_100bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_1one_sub5_100bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_1one_sub5_100bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_1one_sub5_100bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_1one_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub5_100bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffGeneAnn_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub5_100bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub5_100bp,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffCpGann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_1one_sub5_100bp,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffrepeatann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub5_100bp, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub5_100bp, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub5_100bp, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub5_100bp, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_1one_sub5_100bp, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_1one_sub5_100bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_1one_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_1one_sub5_100bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffGeneAnn_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_1one_sub5_100bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub5_100bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffCpGann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_1one_sub5_100bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffrepeatann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub5_100bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub5_100bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub5_100bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub5_100bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_1one_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_1one_sub5_100bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_1one_sub5_100bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_1one_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_1one_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_1one_sub5_100bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_1one_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
#####################################################################################################################













































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
sink()


pdf( file=paste(myOutDir_2two, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two, screeplot=TRUE)
PCASamples(meth_2two)
dev.off()



##################
myDiff_2two = calculateDiffMeth(meth_2two, num.cores=16)
dim(myDiff_2two)
names(myDiff_2two)

myQvalue_2two = myDiff_2two$qvalue
myMethDi_2two = myDiff_2two$meth.diff

pdf(paste(myOutDir_2two,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_2two, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_2two[myQvalue_2two<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_2two,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two[myMethDi_2two>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two[myMethDi_2two>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two[myMethDi_2two<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_2two <- -log10(myQvalue_2two)
myColor_2two[myColor_2two >  2]  = "yes"
myColor_2two[myColor_2two <= 2]  = "no"
DataFrame_Local5_2two <- data.frame(myx1 = -log10(myQvalue_2two), myy1 = myMethDi_2two,   mycolor1 = myColor_2two )

ggplot(DataFrame_Local5_2two, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_2two,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_2two,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_2two  = getMethylDiff(myDiff_2two, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two = getMethylDiff(myDiff_2two, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_2two       = getMethylDiff(myDiff_2two, difference=10, qvalue=0.01)
myDiffTemp_2two      = getMethylDiff(myDiff_2two, difference=0,  qvalue=0.05)

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

write.table(myDiffTemp_2two , file = paste(myOutDir_2two,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
myrepeat.obj = readFeatureFlank(myRepeats,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj = readFeatureFlank(myImprintedRegions1, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj = readFeatureFlank(myImprintedRegions2, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj = readFeatureFlank(myImprintedRegions3, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj = readFeatureFlank(myImprintedRegions4, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj = readFeatureFlank(myImprintedRegions5, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################



########################## annotation for hypo sites.
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
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two,"GRanges"),
                                                 cpg.obj$CpGi,  cpg.obj$shores,
                                                 feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_hyper,  percentage=TRUE)
print(diffCpGann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two,"GRanges"),
                                                    myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                    feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_hyper,  percentage=TRUE)
print(diffrepeatann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two, "GRanges"),
                                                        feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two, "GRanges"),
                                                        feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two, "GRanges"),
                                                        feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two, "GRanges"),
                                                        feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two, "GRanges"),
                                                        feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                        feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_hyper)
sink()

pdf( file=paste(myOutDir_2two, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_2two_sub1_10kb = paste(myOutDir_2two, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_2two_sub1_10kb) ) { dir.create(myOutDir_2two_sub1_10kb, recursive = TRUE) }


tiles_2two_sub1_10kb = tileMethylCounts(filtered.myobj_2two, win.size=10000, step.size=10000)
head(tiles_2two_sub1_10kb[[1]], 3)

meth_2two_sub1_10kb = unite(tiles_2two_sub1_10kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_2two_sub1_10kb)
dim(meth_2two_sub1_10kb)

mat_2two_sub1_10kb = percMethylation(meth_2two_sub1_10kb)
head(mat_2two_sub1_10kb)
dim(mat_2two_sub1_10kb)

write.table(meth_2two_sub1_10kb , 
            file = paste(myOutDir_2two_sub1_10kb,   "0A_meth_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub1_10kb , 
            file = paste(myOutDir_2two_sub1_10kb,   "0B_mat_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two_sub1_10kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two_sub1_10kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub1_10kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two_sub1_10kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two_sub1_10kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub1_10kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_2two_sub1_10kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two_sub1_10kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub1_10kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub1_10kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two_sub1_10kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two_sub1_10kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two_sub1_10kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub1_10kb, screeplot=TRUE)
PCASamples(meth_2two_sub1_10kb)
dev.off()



##################
myDiff_2two_sub1_10kb = calculateDiffMeth(meth_2two_sub1_10kb, num.cores=16)
dim(myDiff_2two_sub1_10kb)
names(myDiff_2two_sub1_10kb)

myQvalue_2two_sub1_10kb = myDiff_2two_sub1_10kb$qvalue
myMethDi_2two_sub1_10kb = myDiff_2two_sub1_10kb$meth.diff

pdf(paste(myOutDir_2two_sub1_10kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_2two_sub1_10kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_2two_sub1_10kb[myQvalue_2two_sub1_10kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_2two_sub1_10kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two_sub1_10kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two_sub1_10kb[myMethDi_2two_sub1_10kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two_sub1_10kb[myMethDi_2two_sub1_10kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two_sub1_10kb[myMethDi_2two_sub1_10kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_2two_sub1_10kb <- -log10(myQvalue_2two_sub1_10kb)
myColor_2two_sub1_10kb[myColor_2two_sub1_10kb >  2]  = "yes"
myColor_2two_sub1_10kb[myColor_2two_sub1_10kb <= 2]  = "no"
DataFrame_Local5_2two_sub1_10kb <- data.frame(myx1 = -log10(myQvalue_2two_sub1_10kb), myy1 = myMethDi_2two_sub1_10kb,   mycolor1 = myColor_2two_sub1_10kb )

ggplot(DataFrame_Local5_2two_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_sub1_10kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_2two_sub1_10kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_2two_sub1_10kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_2two_sub1_10kb  = getMethylDiff(myDiff_2two_sub1_10kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two_sub1_10kb = getMethylDiff(myDiff_2two_sub1_10kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_2two_sub1_10kb       = getMethylDiff(myDiff_2two_sub1_10kb, difference=10, qvalue=0.01)
myDiffTemp_2two_sub1_10kb      = getMethylDiff(myDiff_2two_sub1_10kb, difference=0,  qvalue=0.05)

write.table(myDiff_2two_sub1_10kb , file = paste(myOutDir_2two_sub1_10kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_sub1_10kb , file = paste(myOutDir_2two_sub1_10kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_sub1_10kb , file = paste(myOutDir_2two_sub1_10kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_sub1_10kb , file = paste(myOutDir_2two_sub1_10kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_2two_sub1_10kb , file = paste(myOutDir_2two_sub1_10kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two_sub1_10kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_sub1_10kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two_sub1_10kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_sub1_10kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_2two_sub1_10kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub1_10kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub1_10kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub1_10kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffCpGann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub1_10kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffrepeatann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub1_10kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub1_10kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub1_10kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub1_10kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub1_10kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_sub1_10kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub1_10kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub1_10kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub1_10kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffCpGann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub1_10kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffrepeatann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub1_10kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub1_10kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub1_10kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub1_10kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub1_10kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub1_10kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub1_10kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()







### Tiling windows analysis
##################
myOutDir_2two_sub2_5kb = paste(myOutDir_2two, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_2two_sub2_5kb) ) { dir.create(myOutDir_2two_sub2_5kb, recursive = TRUE) }


tiles_2two_sub2_5kb = tileMethylCounts(filtered.myobj_2two, win.size=5000, step.size=5000)
head(tiles_2two_sub2_5kb[[1]], 3)

meth_2two_sub2_5kb = unite(tiles_2two_sub2_5kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_2two_sub2_5kb)
dim(meth_2two_sub2_5kb)

mat_2two_sub2_5kb = percMethylation(meth_2two_sub2_5kb)
head(mat_2two_sub2_5kb)
dim(mat_2two_sub2_5kb)

write.table(meth_2two_sub2_5kb , 
            file = paste(myOutDir_2two_sub2_5kb,   "0A_meth_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub2_5kb , 
            file = paste(myOutDir_2two_sub2_5kb,   "0B_mat_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two_sub2_5kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two_sub2_5kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub2_5kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two_sub2_5kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two_sub2_5kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub2_5kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_2two_sub2_5kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two_sub2_5kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_5kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_5kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two_sub2_5kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two_sub2_5kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two_sub2_5kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub2_5kb, screeplot=TRUE)
PCASamples(meth_2two_sub2_5kb)
dev.off()



##################
myDiff_2two_sub2_5kb = calculateDiffMeth(meth_2two_sub2_5kb, num.cores=16)
dim(myDiff_2two_sub2_5kb)
names(myDiff_2two_sub2_5kb)

myQvalue_2two_sub2_5kb = myDiff_2two_sub2_5kb$qvalue
myMethDi_2two_sub2_5kb = myDiff_2two_sub2_5kb$meth.diff

pdf(paste(myOutDir_2two_sub2_5kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_2two_sub2_5kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_2two_sub2_5kb[myQvalue_2two_sub2_5kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_2two_sub2_5kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two_sub2_5kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two_sub2_5kb[myMethDi_2two_sub2_5kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two_sub2_5kb[myMethDi_2two_sub2_5kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two_sub2_5kb[myMethDi_2two_sub2_5kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_2two_sub2_5kb <- -log10(myQvalue_2two_sub2_5kb)
myColor_2two_sub2_5kb[myColor_2two_sub2_5kb >  2]  = "yes"
myColor_2two_sub2_5kb[myColor_2two_sub2_5kb <= 2]  = "no"
DataFrame_Local5_2two_sub2_5kb <- data.frame(myx1 = -log10(myQvalue_2two_sub2_5kb), myy1 = myMethDi_2two_sub2_5kb,   mycolor1 = myColor_2two_sub2_5kb )

ggplot(DataFrame_Local5_2two_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_sub2_5kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_2two_sub2_5kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_2two_sub2_5kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_2two_sub2_5kb  = getMethylDiff(myDiff_2two_sub2_5kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two_sub2_5kb = getMethylDiff(myDiff_2two_sub2_5kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_2two_sub2_5kb       = getMethylDiff(myDiff_2two_sub2_5kb, difference=10, qvalue=0.01)
myDiffTemp_2two_sub2_5kb      = getMethylDiff(myDiff_2two_sub2_5kb, difference=0,  qvalue=0.05)

write.table(myDiff_2two_sub2_5kb , file = paste(myOutDir_2two_sub2_5kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_sub2_5kb , file = paste(myOutDir_2two_sub2_5kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_sub2_5kb , file = paste(myOutDir_2two_sub2_5kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_sub2_5kb , file = paste(myOutDir_2two_sub2_5kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_2two_sub2_5kb , file = paste(myOutDir_2two_sub2_5kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two_sub2_5kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_sub2_5kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two_sub2_5kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_sub2_5kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_2two_sub2_5kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub2_5kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub2_5kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2_5kb,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffCpGann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2_5kb,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffrepeatann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_5kb, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_5kb, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_5kb, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_5kb, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_5kb, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_sub2_5kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub2_5kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub2_5kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2_5kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffCpGann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2_5kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffrepeatann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_5kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_5kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_5kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_5kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_5kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub2_5kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub2_5kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()










### Tiling windows analysis
##################
myOutDir_2two_sub3_1kb = paste(myOutDir_2two, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_2two_sub3_1kb) ) { dir.create(myOutDir_2two_sub3_1kb, recursive = TRUE) }


tiles_2two_sub3_1kb = tileMethylCounts(filtered.myobj_2two, win.size=1000, step.size=1000)
head(tiles_2two_sub3_1kb[[1]], 3)

meth_2two_sub3_1kb = unite(tiles_2two_sub3_1kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_2two_sub3_1kb)
dim(meth_2two_sub3_1kb)

mat_2two_sub3_1kb = percMethylation(meth_2two_sub3_1kb)
head(mat_2two_sub3_1kb)
dim(mat_2two_sub3_1kb)

write.table(meth_2two_sub3_1kb , 
            file = paste(myOutDir_2two_sub3_1kb,   "0A_meth_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub3_1kb , 
            file = paste(myOutDir_2two_sub3_1kb,   "0B_mat_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two_sub3_1kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two_sub3_1kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub3_1kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two_sub3_1kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two_sub3_1kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub3_1kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_2two_sub3_1kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two_sub3_1kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub3_1kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub3_1kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two_sub3_1kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two_sub3_1kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two_sub3_1kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub3_1kb, screeplot=TRUE)
PCASamples(meth_2two_sub3_1kb)
dev.off()



##################
myDiff_2two_sub3_1kb = calculateDiffMeth(meth_2two_sub3_1kb, num.cores=16)
dim(myDiff_2two_sub3_1kb)
names(myDiff_2two_sub3_1kb)

myQvalue_2two_sub3_1kb = myDiff_2two_sub3_1kb$qvalue
myMethDi_2two_sub3_1kb = myDiff_2two_sub3_1kb$meth.diff

pdf(paste(myOutDir_2two_sub3_1kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_2two_sub3_1kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_2two_sub3_1kb[myQvalue_2two_sub3_1kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_2two_sub3_1kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two_sub3_1kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two_sub3_1kb[myMethDi_2two_sub3_1kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two_sub3_1kb[myMethDi_2two_sub3_1kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two_sub3_1kb[myMethDi_2two_sub3_1kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_2two_sub3_1kb <- -log10(myQvalue_2two_sub3_1kb)
myColor_2two_sub3_1kb[myColor_2two_sub3_1kb >  2]  = "yes"
myColor_2two_sub3_1kb[myColor_2two_sub3_1kb <= 2]  = "no"
DataFrame_Local5_2two_sub3_1kb <- data.frame(myx1 = -log10(myQvalue_2two_sub3_1kb), myy1 = myMethDi_2two_sub3_1kb,   mycolor1 = myColor_2two_sub3_1kb )

ggplot(DataFrame_Local5_2two_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_sub3_1kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_2two_sub3_1kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_2two_sub3_1kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_2two_sub3_1kb  = getMethylDiff(myDiff_2two_sub3_1kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two_sub3_1kb = getMethylDiff(myDiff_2two_sub3_1kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_2two_sub3_1kb       = getMethylDiff(myDiff_2two_sub3_1kb, difference=10, qvalue=0.01)
myDiffTemp_2two_sub3_1kb      = getMethylDiff(myDiff_2two_sub3_1kb, difference=0,  qvalue=0.05)

write.table(myDiff_2two_sub3_1kb , file = paste(myOutDir_2two_sub3_1kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_sub3_1kb , file = paste(myOutDir_2two_sub3_1kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_sub3_1kb , file = paste(myOutDir_2two_sub3_1kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_sub3_1kb , file = paste(myOutDir_2two_sub3_1kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_2two_sub3_1kb , file = paste(myOutDir_2two_sub3_1kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two_sub3_1kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_sub3_1kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two_sub3_1kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_sub3_1kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_2two_sub3_1kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub3_1kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub3_1kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub3_1kb,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffCpGann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub3_1kb,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffrepeatann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub3_1kb, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub3_1kb, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub3_1kb, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub3_1kb, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub3_1kb, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_sub3_1kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub3_1kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub3_1kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub3_1kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffCpGann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub3_1kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffrepeatann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub3_1kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub3_1kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub3_1kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub3_1kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub3_1kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub3_1kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub3_1kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_2two_sub4_500bp = paste(myOutDir_2two, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_2two_sub4_500bp) ) { dir.create(myOutDir_2two_sub4_500bp, recursive = TRUE) }


tiles_2two_sub4_500bp = tileMethylCounts(filtered.myobj_2two, win.size=500, step.size=500)
head(tiles_2two_sub4_500bp[[1]], 3)

meth_2two_sub4_500bp = unite(tiles_2two_sub4_500bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_2two_sub4_500bp)
dim(meth_2two_sub4_500bp)

mat_2two_sub4_500bp = percMethylation(meth_2two_sub4_500bp)
head(mat_2two_sub4_500bp)
dim(mat_2two_sub4_500bp)

write.table(meth_2two_sub4_500bp , 
            file = paste(myOutDir_2two_sub4_500bp,   "0A_meth_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub4_500bp , 
            file = paste(myOutDir_2two_sub4_500bp,   "0B_mat_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two_sub4_500bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two_sub4_500bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub4_500bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two_sub4_500bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two_sub4_500bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub4_500bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_2two_sub4_500bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two_sub4_500bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub4_500bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub4_500bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two_sub4_500bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two_sub4_500bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two_sub4_500bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub4_500bp, screeplot=TRUE)
PCASamples(meth_2two_sub4_500bp)
dev.off()



##################
myDiff_2two_sub4_500bp = calculateDiffMeth(meth_2two_sub4_500bp, num.cores=16)
dim(myDiff_2two_sub4_500bp)
names(myDiff_2two_sub4_500bp)

myQvalue_2two_sub4_500bp = myDiff_2two_sub4_500bp$qvalue
myMethDi_2two_sub4_500bp = myDiff_2two_sub4_500bp$meth.diff

pdf(paste(myOutDir_2two_sub4_500bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_2two_sub4_500bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_2two_sub4_500bp[myQvalue_2two_sub4_500bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_2two_sub4_500bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two_sub4_500bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two_sub4_500bp[myMethDi_2two_sub4_500bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two_sub4_500bp[myMethDi_2two_sub4_500bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two_sub4_500bp[myMethDi_2two_sub4_500bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_2two_sub4_500bp <- -log10(myQvalue_2two_sub4_500bp)
myColor_2two_sub4_500bp[myColor_2two_sub4_500bp >  2]  = "yes"
myColor_2two_sub4_500bp[myColor_2two_sub4_500bp <= 2]  = "no"
DataFrame_Local5_2two_sub4_500bp <- data.frame(myx1 = -log10(myQvalue_2two_sub4_500bp), myy1 = myMethDi_2two_sub4_500bp,   mycolor1 = myColor_2two_sub4_500bp )

ggplot(DataFrame_Local5_2two_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_sub4_500bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_2two_sub4_500bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_2two_sub4_500bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_2two_sub4_500bp  = getMethylDiff(myDiff_2two_sub4_500bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two_sub4_500bp = getMethylDiff(myDiff_2two_sub4_500bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_2two_sub4_500bp       = getMethylDiff(myDiff_2two_sub4_500bp, difference=10, qvalue=0.01)
myDiffTemp_2two_sub4_500bp      = getMethylDiff(myDiff_2two_sub4_500bp, difference=0,  qvalue=0.05)

write.table(myDiff_2two_sub4_500bp , file = paste(myOutDir_2two_sub4_500bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_sub4_500bp , file = paste(myOutDir_2two_sub4_500bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_sub4_500bp , file = paste(myOutDir_2two_sub4_500bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_sub4_500bp , file = paste(myOutDir_2two_sub4_500bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_2two_sub4_500bp , file = paste(myOutDir_2two_sub4_500bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two_sub4_500bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_sub4_500bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two_sub4_500bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_sub4_500bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_2two_sub4_500bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub4_500bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub4_500bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub4_500bp,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffCpGann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub4_500bp,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffrepeatann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub4_500bp, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub4_500bp, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub4_500bp, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub4_500bp, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub4_500bp, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_sub4_500bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub4_500bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub4_500bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub4_500bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffCpGann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub4_500bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffrepeatann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub4_500bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub4_500bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub4_500bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub4_500bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub4_500bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub4_500bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub4_500bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
















### Tiling windows analysis
##################
myOutDir_2two_sub5_100bp = paste(myOutDir_2two, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_2two_sub5_100bp) ) { dir.create(myOutDir_2two_sub5_100bp, recursive = TRUE) }


tiles_2two_sub5_100bp = tileMethylCounts(filtered.myobj_2two, win.size=100, step.size=100)
head(tiles_2two_sub5_100bp[[1]], 3)

meth_2two_sub5_100bp = unite(tiles_2two_sub5_100bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_2two_sub5_100bp)
dim(meth_2two_sub5_100bp)

mat_2two_sub5_100bp = percMethylation(meth_2two_sub5_100bp)
head(mat_2two_sub5_100bp)
dim(mat_2two_sub5_100bp)

write.table(meth_2two_sub5_100bp , 
            file = paste(myOutDir_2two_sub5_100bp,   "0A_meth_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub5_100bp , 
            file = paste(myOutDir_2two_sub5_100bp,   "0B_mat_2two_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two_sub5_100bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two_sub5_100bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub5_100bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two_sub5_100bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two_sub5_100bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub5_100bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_2two_sub5_100bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_2two_sub5_100bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub5_100bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub5_100bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_2two_sub5_100bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_2two_sub5_100bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_2two_sub5_100bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_2two_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub5_100bp, screeplot=TRUE)
PCASamples(meth_2two_sub5_100bp)
dev.off()



##################
myDiff_2two_sub5_100bp = calculateDiffMeth(meth_2two_sub5_100bp, num.cores=16)
dim(myDiff_2two_sub5_100bp)
names(myDiff_2two_sub5_100bp)

myQvalue_2two_sub5_100bp = myDiff_2two_sub5_100bp$qvalue
myMethDi_2two_sub5_100bp = myDiff_2two_sub5_100bp$meth.diff

pdf(paste(myOutDir_2two_sub5_100bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_2two_sub5_100bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_2two_sub5_100bp[myQvalue_2two_sub5_100bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_2two_sub5_100bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two_sub5_100bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two_sub5_100bp[myMethDi_2two_sub5_100bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two_sub5_100bp[myMethDi_2two_sub5_100bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two_sub5_100bp[myMethDi_2two_sub5_100bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_2two_sub5_100bp <- -log10(myQvalue_2two_sub5_100bp)
myColor_2two_sub5_100bp[myColor_2two_sub5_100bp >  2]  = "yes"
myColor_2two_sub5_100bp[myColor_2two_sub5_100bp <= 2]  = "no"
DataFrame_Local5_2two_sub5_100bp <- data.frame(myx1 = -log10(myQvalue_2two_sub5_100bp), myy1 = myMethDi_2two_sub5_100bp,   mycolor1 = myColor_2two_sub5_100bp )

ggplot(DataFrame_Local5_2two_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_sub5_100bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_2two_sub5_100bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_2two_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_2two_sub5_100bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_2two_sub5_100bp  = getMethylDiff(myDiff_2two_sub5_100bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two_sub5_100bp = getMethylDiff(myDiff_2two_sub5_100bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_2two_sub5_100bp       = getMethylDiff(myDiff_2two_sub5_100bp, difference=10, qvalue=0.01)
myDiffTemp_2two_sub5_100bp      = getMethylDiff(myDiff_2two_sub5_100bp, difference=0,  qvalue=0.05)

write.table(myDiff_2two_sub5_100bp , file = paste(myOutDir_2two_sub5_100bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_sub5_100bp , file = paste(myOutDir_2two_sub5_100bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_sub5_100bp , file = paste(myOutDir_2two_sub5_100bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_sub5_100bp , file = paste(myOutDir_2two_sub5_100bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_2two_sub5_100bp , file = paste(myOutDir_2two_sub5_100bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two_sub5_100bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_sub5_100bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two_sub5_100bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_sub5_100bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_2two_sub5_100bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub5_100bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub5_100bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub5_100bp,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffCpGann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub5_100bp,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffrepeatann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub5_100bp, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub5_100bp, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub5_100bp, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub5_100bp, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub5_100bp, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_sub5_100bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_sub5_100bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub5_100bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub5_100bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffCpGann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub5_100bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffrepeatann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub5_100bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub5_100bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub5_100bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub5_100bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub5_100bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_sub5_100bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_2two_sub5_100bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
#####################################################################################################################



























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
sink()


pdf( file=paste(myOutDir_3three, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three, screeplot=TRUE)
PCASamples(meth_3three)
dev.off()



##################
myDiff_3three = calculateDiffMeth(meth_3three, num.cores=16)
dim(myDiff_3three)
names(myDiff_3three)

myQvalue_3three = myDiff_3three$qvalue
myMethDi_3three = myDiff_3three$meth.diff

pdf(paste(myOutDir_3three,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_3three, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_3three[myQvalue_3three<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_3three,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three[myMethDi_3three>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three[myMethDi_3three>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three[myMethDi_3three<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_3three <- -log10(myQvalue_3three)
myColor_3three[myColor_3three >  2]  = "yes"
myColor_3three[myColor_3three <= 2]  = "no"
DataFrame_Local5_3three <- data.frame(myx1 = -log10(myQvalue_3three), myy1 = myMethDi_3three,   mycolor1 = myColor_3three )

ggplot(DataFrame_Local5_3three, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_3three,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_3three,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_3three  = getMethylDiff(myDiff_3three, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three = getMethylDiff(myDiff_3three, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_3three       = getMethylDiff(myDiff_3three, difference=10, qvalue=0.01)
myDiffTemp_3three      = getMethylDiff(myDiff_3three, difference=0,  qvalue=0.05)

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

write.table(myDiffTemp_3three , file = paste(myOutDir_3three,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
myrepeat.obj = readFeatureFlank(myRepeats,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj = readFeatureFlank(myImprintedRegions1, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj = readFeatureFlank(myImprintedRegions2, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj = readFeatureFlank(myImprintedRegions3, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj = readFeatureFlank(myImprintedRegions4, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj = readFeatureFlank(myImprintedRegions5, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################



########################## annotation for hypo sites.
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
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three,"GRanges"),
                                                   cpg.obj$CpGi,  cpg.obj$shores,
                                                   feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_hyper,  percentage=TRUE)
print(diffCpGann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three,"GRanges"),
                                                      myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                      feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_hyper,  percentage=TRUE)
print(diffrepeatann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three, "GRanges"),
                                                          feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                          feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three, "GRanges"),
                                                          feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                          feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three, "GRanges"),
                                                          feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                          feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three, "GRanges"),
                                                          feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                          feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three, "GRanges"),
                                                          feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                          feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_hyper)
sink()

pdf( file=paste(myOutDir_3three, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_3three_sub1_10kb = paste(myOutDir_3three, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_3three_sub1_10kb) ) { dir.create(myOutDir_3three_sub1_10kb, recursive = TRUE) }


tiles_3three_sub1_10kb = tileMethylCounts(filtered.myobj_3three, win.size=10000, step.size=10000)
head(tiles_3three_sub1_10kb[[1]], 3)

meth_3three_sub1_10kb = unite(tiles_3three_sub1_10kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_3three_sub1_10kb)
dim(meth_3three_sub1_10kb)

mat_3three_sub1_10kb = percMethylation(meth_3three_sub1_10kb)
head(mat_3three_sub1_10kb)
dim(mat_3three_sub1_10kb)

write.table(meth_3three_sub1_10kb , 
            file = paste(myOutDir_3three_sub1_10kb,   "0A_meth_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub1_10kb , 
            file = paste(myOutDir_3three_sub1_10kb,   "0B_mat_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three_sub1_10kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three_sub1_10kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub1_10kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three_sub1_10kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three_sub1_10kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub1_10kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_3three_sub1_10kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three_sub1_10kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub1_10kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub1_10kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three_sub1_10kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three_sub1_10kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three_sub1_10kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub1_10kb, screeplot=TRUE)
PCASamples(meth_3three_sub1_10kb)
dev.off()



##################
myDiff_3three_sub1_10kb = calculateDiffMeth(meth_3three_sub1_10kb, num.cores=16)
dim(myDiff_3three_sub1_10kb)
names(myDiff_3three_sub1_10kb)

myQvalue_3three_sub1_10kb = myDiff_3three_sub1_10kb$qvalue
myMethDi_3three_sub1_10kb = myDiff_3three_sub1_10kb$meth.diff

pdf(paste(myOutDir_3three_sub1_10kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_3three_sub1_10kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_3three_sub1_10kb[myQvalue_3three_sub1_10kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_3three_sub1_10kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three_sub1_10kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three_sub1_10kb[myMethDi_3three_sub1_10kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three_sub1_10kb[myMethDi_3three_sub1_10kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three_sub1_10kb[myMethDi_3three_sub1_10kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_3three_sub1_10kb <- -log10(myQvalue_3three_sub1_10kb)
myColor_3three_sub1_10kb[myColor_3three_sub1_10kb >  2]  = "yes"
myColor_3three_sub1_10kb[myColor_3three_sub1_10kb <= 2]  = "no"
DataFrame_Local5_3three_sub1_10kb <- data.frame(myx1 = -log10(myQvalue_3three_sub1_10kb), myy1 = myMethDi_3three_sub1_10kb,   mycolor1 = myColor_3three_sub1_10kb )

ggplot(DataFrame_Local5_3three_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_sub1_10kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_3three_sub1_10kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_3three_sub1_10kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_3three_sub1_10kb  = getMethylDiff(myDiff_3three_sub1_10kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three_sub1_10kb = getMethylDiff(myDiff_3three_sub1_10kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_3three_sub1_10kb       = getMethylDiff(myDiff_3three_sub1_10kb, difference=10, qvalue=0.01)
myDiffTemp_3three_sub1_10kb      = getMethylDiff(myDiff_3three_sub1_10kb, difference=0,  qvalue=0.05)

write.table(myDiff_3three_sub1_10kb , file = paste(myOutDir_3three_sub1_10kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_sub1_10kb , file = paste(myOutDir_3three_sub1_10kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_sub1_10kb , file = paste(myOutDir_3three_sub1_10kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_sub1_10kb , file = paste(myOutDir_3three_sub1_10kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_3three_sub1_10kb , file = paste(myOutDir_3three_sub1_10kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three_sub1_10kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_sub1_10kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three_sub1_10kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_sub1_10kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_3three_sub1_10kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub1_10kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub1_10kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub1_10kb,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffCpGann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub1_10kb,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffrepeatann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub1_10kb, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub1_10kb, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub1_10kb, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub1_10kb, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub1_10kb, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_sub1_10kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub1_10kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub1_10kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub1_10kb,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffCpGann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub1_10kb,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffrepeatann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub1_10kb, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub1_10kb, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub1_10kb, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub1_10kb, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub1_10kb, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub1_10kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub1_10kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()







### Tiling windows analysis
##################
myOutDir_3three_sub2_5kb = paste(myOutDir_3three, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_3three_sub2_5kb) ) { dir.create(myOutDir_3three_sub2_5kb, recursive = TRUE) }


tiles_3three_sub2_5kb = tileMethylCounts(filtered.myobj_3three, win.size=5000, step.size=5000)
head(tiles_3three_sub2_5kb[[1]], 3)

meth_3three_sub2_5kb = unite(tiles_3three_sub2_5kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_3three_sub2_5kb)
dim(meth_3three_sub2_5kb)

mat_3three_sub2_5kb = percMethylation(meth_3three_sub2_5kb)
head(mat_3three_sub2_5kb)
dim(mat_3three_sub2_5kb)

write.table(meth_3three_sub2_5kb , 
            file = paste(myOutDir_3three_sub2_5kb,   "0A_meth_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub2_5kb , 
            file = paste(myOutDir_3three_sub2_5kb,   "0B_mat_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three_sub2_5kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three_sub2_5kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub2_5kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three_sub2_5kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three_sub2_5kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub2_5kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_3three_sub2_5kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three_sub2_5kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub2_5kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub2_5kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three_sub2_5kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three_sub2_5kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three_sub2_5kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub2_5kb, screeplot=TRUE)
PCASamples(meth_3three_sub2_5kb)
dev.off()



##################
myDiff_3three_sub2_5kb = calculateDiffMeth(meth_3three_sub2_5kb, num.cores=16)
dim(myDiff_3three_sub2_5kb)
names(myDiff_3three_sub2_5kb)

myQvalue_3three_sub2_5kb = myDiff_3three_sub2_5kb$qvalue
myMethDi_3three_sub2_5kb = myDiff_3three_sub2_5kb$meth.diff

pdf(paste(myOutDir_3three_sub2_5kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_3three_sub2_5kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_3three_sub2_5kb[myQvalue_3three_sub2_5kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_3three_sub2_5kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three_sub2_5kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three_sub2_5kb[myMethDi_3three_sub2_5kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three_sub2_5kb[myMethDi_3three_sub2_5kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three_sub2_5kb[myMethDi_3three_sub2_5kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_3three_sub2_5kb <- -log10(myQvalue_3three_sub2_5kb)
myColor_3three_sub2_5kb[myColor_3three_sub2_5kb >  2]  = "yes"
myColor_3three_sub2_5kb[myColor_3three_sub2_5kb <= 2]  = "no"
DataFrame_Local5_3three_sub2_5kb <- data.frame(myx1 = -log10(myQvalue_3three_sub2_5kb), myy1 = myMethDi_3three_sub2_5kb,   mycolor1 = myColor_3three_sub2_5kb )

ggplot(DataFrame_Local5_3three_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_sub2_5kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_3three_sub2_5kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_3three_sub2_5kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_3three_sub2_5kb  = getMethylDiff(myDiff_3three_sub2_5kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three_sub2_5kb = getMethylDiff(myDiff_3three_sub2_5kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_3three_sub2_5kb       = getMethylDiff(myDiff_3three_sub2_5kb, difference=10, qvalue=0.01)
myDiffTemp_3three_sub2_5kb      = getMethylDiff(myDiff_3three_sub2_5kb, difference=0,  qvalue=0.05)

write.table(myDiff_3three_sub2_5kb , file = paste(myOutDir_3three_sub2_5kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_sub2_5kb , file = paste(myOutDir_3three_sub2_5kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_sub2_5kb , file = paste(myOutDir_3three_sub2_5kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_sub2_5kb , file = paste(myOutDir_3three_sub2_5kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_3three_sub2_5kb , file = paste(myOutDir_3three_sub2_5kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three_sub2_5kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_sub2_5kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three_sub2_5kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_sub2_5kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_3three_sub2_5kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub2_5kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub2_5kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub2_5kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffCpGann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub2_5kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffrepeatann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub2_5kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub2_5kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub2_5kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub2_5kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub2_5kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_sub2_5kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub2_5kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub2_5kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub2_5kb,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffCpGann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub2_5kb,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffrepeatann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub2_5kb, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub2_5kb, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub2_5kb, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub2_5kb, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub2_5kb, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub2_5kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub2_5kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()










### Tiling windows analysis
##################
myOutDir_3three_sub3_1kb = paste(myOutDir_3three, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_3three_sub3_1kb) ) { dir.create(myOutDir_3three_sub3_1kb, recursive = TRUE) }


tiles_3three_sub3_1kb = tileMethylCounts(filtered.myobj_3three, win.size=1000, step.size=1000)
head(tiles_3three_sub3_1kb[[1]], 3)

meth_3three_sub3_1kb = unite(tiles_3three_sub3_1kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_3three_sub3_1kb)
dim(meth_3three_sub3_1kb)

mat_3three_sub3_1kb = percMethylation(meth_3three_sub3_1kb)
head(mat_3three_sub3_1kb)
dim(mat_3three_sub3_1kb)

write.table(meth_3three_sub3_1kb , 
            file = paste(myOutDir_3three_sub3_1kb,   "0A_meth_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub3_1kb , 
            file = paste(myOutDir_3three_sub3_1kb,   "0B_mat_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three_sub3_1kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three_sub3_1kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub3_1kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three_sub3_1kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three_sub3_1kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub3_1kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_3three_sub3_1kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three_sub3_1kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_1kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_1kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three_sub3_1kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three_sub3_1kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three_sub3_1kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub3_1kb, screeplot=TRUE)
PCASamples(meth_3three_sub3_1kb)
dev.off()



##################
myDiff_3three_sub3_1kb = calculateDiffMeth(meth_3three_sub3_1kb, num.cores=16)
dim(myDiff_3three_sub3_1kb)
names(myDiff_3three_sub3_1kb)

myQvalue_3three_sub3_1kb = myDiff_3three_sub3_1kb$qvalue
myMethDi_3three_sub3_1kb = myDiff_3three_sub3_1kb$meth.diff

pdf(paste(myOutDir_3three_sub3_1kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_3three_sub3_1kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_3three_sub3_1kb[myQvalue_3three_sub3_1kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_3three_sub3_1kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three_sub3_1kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_3three_sub3_1kb <- -log10(myQvalue_3three_sub3_1kb)
myColor_3three_sub3_1kb[myColor_3three_sub3_1kb >  2]  = "yes"
myColor_3three_sub3_1kb[myColor_3three_sub3_1kb <= 2]  = "no"
DataFrame_Local5_3three_sub3_1kb <- data.frame(myx1 = -log10(myQvalue_3three_sub3_1kb), myy1 = myMethDi_3three_sub3_1kb,   mycolor1 = myColor_3three_sub3_1kb )

ggplot(DataFrame_Local5_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_sub3_1kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_3three_sub3_1kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_3three_sub3_1kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_3three_sub3_1kb  = getMethylDiff(myDiff_3three_sub3_1kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three_sub3_1kb = getMethylDiff(myDiff_3three_sub3_1kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_3three_sub3_1kb       = getMethylDiff(myDiff_3three_sub3_1kb, difference=10, qvalue=0.01)
myDiffTemp_3three_sub3_1kb      = getMethylDiff(myDiff_3three_sub3_1kb, difference=0,  qvalue=0.05)

write.table(myDiff_3three_sub3_1kb , file = paste(myOutDir_3three_sub3_1kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_sub3_1kb , file = paste(myOutDir_3three_sub3_1kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_sub3_1kb , file = paste(myOutDir_3three_sub3_1kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_sub3_1kb , file = paste(myOutDir_3three_sub3_1kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_3three_sub3_1kb , file = paste(myOutDir_3three_sub3_1kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three_sub3_1kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_sub3_1kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three_sub3_1kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_sub3_1kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_3three_sub3_1kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub3_1kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub3_1kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub3_1kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffCpGann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub3_1kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffrepeatann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_sub3_1kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub3_1kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub3_1kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub3_1kb,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffCpGann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub3_1kb,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffrepeatann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub3_1kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub3_1kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_3three_sub4_500bp = paste(myOutDir_3three, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_3three_sub4_500bp) ) { dir.create(myOutDir_3three_sub4_500bp, recursive = TRUE) }


tiles_3three_sub4_500bp = tileMethylCounts(filtered.myobj_3three, win.size=500, step.size=500)
head(tiles_3three_sub4_500bp[[1]], 3)

meth_3three_sub4_500bp = unite(tiles_3three_sub4_500bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_3three_sub4_500bp)
dim(meth_3three_sub4_500bp)

mat_3three_sub4_500bp = percMethylation(meth_3three_sub4_500bp)
head(mat_3three_sub4_500bp)
dim(mat_3three_sub4_500bp)

write.table(meth_3three_sub4_500bp , 
            file = paste(myOutDir_3three_sub4_500bp,   "0A_meth_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub4_500bp , 
            file = paste(myOutDir_3three_sub4_500bp,   "0B_mat_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three_sub4_500bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three_sub4_500bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub4_500bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three_sub4_500bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three_sub4_500bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub4_500bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_3three_sub4_500bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three_sub4_500bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub4_500bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub4_500bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three_sub4_500bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three_sub4_500bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three_sub4_500bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub4_500bp, screeplot=TRUE)
PCASamples(meth_3three_sub4_500bp)
dev.off()



##################
myDiff_3three_sub4_500bp = calculateDiffMeth(meth_3three_sub4_500bp, num.cores=16)
dim(myDiff_3three_sub4_500bp)
names(myDiff_3three_sub4_500bp)

myQvalue_3three_sub4_500bp = myDiff_3three_sub4_500bp$qvalue
myMethDi_3three_sub4_500bp = myDiff_3three_sub4_500bp$meth.diff

pdf(paste(myOutDir_3three_sub4_500bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_3three_sub4_500bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_3three_sub4_500bp[myQvalue_3three_sub4_500bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_3three_sub4_500bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three_sub4_500bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three_sub4_500bp[myMethDi_3three_sub4_500bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three_sub4_500bp[myMethDi_3three_sub4_500bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three_sub4_500bp[myMethDi_3three_sub4_500bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_3three_sub4_500bp <- -log10(myQvalue_3three_sub4_500bp)
myColor_3three_sub4_500bp[myColor_3three_sub4_500bp >  2]  = "yes"
myColor_3three_sub4_500bp[myColor_3three_sub4_500bp <= 2]  = "no"
DataFrame_Local5_3three_sub4_500bp <- data.frame(myx1 = -log10(myQvalue_3three_sub4_500bp), myy1 = myMethDi_3three_sub4_500bp,   mycolor1 = myColor_3three_sub4_500bp )

ggplot(DataFrame_Local5_3three_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_sub4_500bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_3three_sub4_500bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_3three_sub4_500bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_3three_sub4_500bp  = getMethylDiff(myDiff_3three_sub4_500bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three_sub4_500bp = getMethylDiff(myDiff_3three_sub4_500bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_3three_sub4_500bp       = getMethylDiff(myDiff_3three_sub4_500bp, difference=10, qvalue=0.01)
myDiffTemp_3three_sub4_500bp      = getMethylDiff(myDiff_3three_sub4_500bp, difference=0,  qvalue=0.05)

write.table(myDiff_3three_sub4_500bp , file = paste(myOutDir_3three_sub4_500bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_sub4_500bp , file = paste(myOutDir_3three_sub4_500bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_sub4_500bp , file = paste(myOutDir_3three_sub4_500bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_sub4_500bp , file = paste(myOutDir_3three_sub4_500bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_3three_sub4_500bp , file = paste(myOutDir_3three_sub4_500bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three_sub4_500bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_sub4_500bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three_sub4_500bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_sub4_500bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_3three_sub4_500bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub4_500bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub4_500bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub4_500bp,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffCpGann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub4_500bp,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffrepeatann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub4_500bp, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub4_500bp, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub4_500bp, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub4_500bp, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub4_500bp, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_sub4_500bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub4_500bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub4_500bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub4_500bp,"GRanges"),
                                                              cpg.obj$CpGi,  cpg.obj$shores,
                                                              feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffCpGann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub4_500bp,"GRanges"),
                                                                 myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                 feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffrepeatann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub4_500bp, "GRanges"),
                                                                     feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub4_500bp, "GRanges"),
                                                                     feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub4_500bp, "GRanges"),
                                                                     feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub4_500bp, "GRanges"),
                                                                     feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub4_500bp, "GRanges"),
                                                                     feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub4_500bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub4_500bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
















### Tiling windows analysis
##################
myOutDir_3three_sub5_100bp = paste(myOutDir_3three, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_3three_sub5_100bp) ) { dir.create(myOutDir_3three_sub5_100bp, recursive = TRUE) }


tiles_3three_sub5_100bp = tileMethylCounts(filtered.myobj_3three, win.size=100, step.size=100)
head(tiles_3three_sub5_100bp[[1]], 3)

meth_3three_sub5_100bp = unite(tiles_3three_sub5_100bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_3three_sub5_100bp)
dim(meth_3three_sub5_100bp)

mat_3three_sub5_100bp = percMethylation(meth_3three_sub5_100bp)
head(mat_3three_sub5_100bp)
dim(mat_3three_sub5_100bp)

write.table(meth_3three_sub5_100bp , 
            file = paste(myOutDir_3three_sub5_100bp,   "0A_meth_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub5_100bp , 
            file = paste(myOutDir_3three_sub5_100bp,   "0B_mat_3three_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three_sub5_100bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three_sub5_100bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub5_100bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three_sub5_100bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three_sub5_100bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub5_100bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_3three_sub5_100bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_3three_sub5_100bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub5_100bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub5_100bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_3three_sub5_100bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_3three_sub5_100bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_3three_sub5_100bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_3three_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub5_100bp, screeplot=TRUE)
PCASamples(meth_3three_sub5_100bp)
dev.off()



##################
myDiff_3three_sub5_100bp = calculateDiffMeth(meth_3three_sub5_100bp, num.cores=16)
dim(myDiff_3three_sub5_100bp)
names(myDiff_3three_sub5_100bp)

myQvalue_3three_sub5_100bp = myDiff_3three_sub5_100bp$qvalue
myMethDi_3three_sub5_100bp = myDiff_3three_sub5_100bp$meth.diff

pdf(paste(myOutDir_3three_sub5_100bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_3three_sub5_100bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_3three_sub5_100bp[myQvalue_3three_sub5_100bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_3three_sub5_100bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three_sub5_100bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three_sub5_100bp[myMethDi_3three_sub5_100bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three_sub5_100bp[myMethDi_3three_sub5_100bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three_sub5_100bp[myMethDi_3three_sub5_100bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_3three_sub5_100bp <- -log10(myQvalue_3three_sub5_100bp)
myColor_3three_sub5_100bp[myColor_3three_sub5_100bp >  2]  = "yes"
myColor_3three_sub5_100bp[myColor_3three_sub5_100bp <= 2]  = "no"
DataFrame_Local5_3three_sub5_100bp <- data.frame(myx1 = -log10(myQvalue_3three_sub5_100bp), myy1 = myMethDi_3three_sub5_100bp,   mycolor1 = myColor_3three_sub5_100bp )

ggplot(DataFrame_Local5_3three_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_sub5_100bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_3three_sub5_100bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_3three_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_3three_sub5_100bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_3three_sub5_100bp  = getMethylDiff(myDiff_3three_sub5_100bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three_sub5_100bp = getMethylDiff(myDiff_3three_sub5_100bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_3three_sub5_100bp       = getMethylDiff(myDiff_3three_sub5_100bp, difference=10, qvalue=0.01)
myDiffTemp_3three_sub5_100bp      = getMethylDiff(myDiff_3three_sub5_100bp, difference=0,  qvalue=0.05)

write.table(myDiff_3three_sub5_100bp , file = paste(myOutDir_3three_sub5_100bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_sub5_100bp , file = paste(myOutDir_3three_sub5_100bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_sub5_100bp , file = paste(myOutDir_3three_sub5_100bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_sub5_100bp , file = paste(myOutDir_3three_sub5_100bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_3three_sub5_100bp , file = paste(myOutDir_3three_sub5_100bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three_sub5_100bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_sub5_100bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three_sub5_100bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_sub5_100bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_3three_sub5_100bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub5_100bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub5_100bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub5_100bp,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffCpGann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub5_100bp,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffrepeatann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub5_100bp, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub5_100bp, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub5_100bp, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub5_100bp, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub5_100bp, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_sub5_100bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_sub5_100bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub5_100bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub5_100bp,"GRanges"),
                                                              cpg.obj$CpGi,  cpg.obj$shores,
                                                              feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffCpGann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub5_100bp,"GRanges"),
                                                                 myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                 feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffrepeatann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub5_100bp, "GRanges"),
                                                                     feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub5_100bp, "GRanges"),
                                                                     feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub5_100bp, "GRanges"),
                                                                     feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub5_100bp, "GRanges"),
                                                                     feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub5_100bp, "GRanges"),
                                                                     feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                     feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_sub5_100bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_3three_sub5_100bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
#####################################################################################################################



































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
sink()


pdf( file=paste(myOutDir_4four, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four, screeplot=TRUE)
PCASamples(meth_4four)
dev.off()



##################
myDiff_4four = calculateDiffMeth(meth_4four, num.cores=16)
dim(myDiff_4four)
names(myDiff_4four)

myQvalue_4four = myDiff_4four$qvalue
myMethDi_4four = myDiff_4four$meth.diff

pdf(paste(myOutDir_4four,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_4four, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_4four[myQvalue_4four<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_4four,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four[myMethDi_4four>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four[myMethDi_4four>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four[myMethDi_4four<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_4four <- -log10(myQvalue_4four)
myColor_4four[myColor_4four >  2]  = "yes"
myColor_4four[myColor_4four <= 2]  = "no"
DataFrame_Local5_4four <- data.frame(myx1 = -log10(myQvalue_4four), myy1 = myMethDi_4four,   mycolor1 = myColor_4four )

ggplot(DataFrame_Local5_4four, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_4four,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_4four,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_4four  = getMethylDiff(myDiff_4four, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four = getMethylDiff(myDiff_4four, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_4four       = getMethylDiff(myDiff_4four, difference=10, qvalue=0.01)
myDiffTemp_4four      = getMethylDiff(myDiff_4four, difference=0,  qvalue=0.05)

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

write.table(myDiffTemp_4four , file = paste(myOutDir_4four,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
myrepeat.obj = readFeatureFlank(myRepeats,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj = readFeatureFlank(myImprintedRegions1, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj = readFeatureFlank(myImprintedRegions2, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj = readFeatureFlank(myImprintedRegions3, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj = readFeatureFlank(myImprintedRegions4, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj = readFeatureFlank(myImprintedRegions5, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################



########################## annotation for hypo sites.
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
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four,"GRanges"),
                                                  cpg.obj$CpGi,  cpg.obj$shores,
                                                  feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_hyper,  percentage=TRUE)
print(diffCpGann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four,"GRanges"),
                                                     myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                     feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_hyper,  percentage=TRUE)
print(diffrepeatann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four, "GRanges"),
                                                         feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four, "GRanges"),
                                                         feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four, "GRanges"),
                                                         feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four, "GRanges"),
                                                         feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four, "GRanges"),
                                                         feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                         feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_hyper)
sink()

pdf( file=paste(myOutDir_4four, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_4four_sub1_10kb = paste(myOutDir_4four, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_4four_sub1_10kb) ) { dir.create(myOutDir_4four_sub1_10kb, recursive = TRUE) }


tiles_4four_sub1_10kb = tileMethylCounts(filtered.myobj_4four, win.size=10000, step.size=10000)
head(tiles_4four_sub1_10kb[[1]], 3)

meth_4four_sub1_10kb = unite(tiles_4four_sub1_10kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_4four_sub1_10kb)
dim(meth_4four_sub1_10kb)

mat_4four_sub1_10kb = percMethylation(meth_4four_sub1_10kb)
head(mat_4four_sub1_10kb)
dim(mat_4four_sub1_10kb)

write.table(meth_4four_sub1_10kb , 
            file = paste(myOutDir_4four_sub1_10kb,   "0A_meth_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub1_10kb , 
            file = paste(myOutDir_4four_sub1_10kb,   "0B_mat_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four_sub1_10kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four_sub1_10kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub1_10kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four_sub1_10kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub1_10kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four_sub1_10kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub1_10kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_4four_sub1_10kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four_sub1_10kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub1_10kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub1_10kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four_sub1_10kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four_sub1_10kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four_sub1_10kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub1_10kb, screeplot=TRUE)
PCASamples(meth_4four_sub1_10kb)
dev.off()



##################
myDiff_4four_sub1_10kb = calculateDiffMeth(meth_4four_sub1_10kb, num.cores=16)
dim(myDiff_4four_sub1_10kb)
names(myDiff_4four_sub1_10kb)

myQvalue_4four_sub1_10kb = myDiff_4four_sub1_10kb$qvalue
myMethDi_4four_sub1_10kb = myDiff_4four_sub1_10kb$meth.diff

pdf(paste(myOutDir_4four_sub1_10kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_4four_sub1_10kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_4four_sub1_10kb[myQvalue_4four_sub1_10kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_4four_sub1_10kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four_sub1_10kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four_sub1_10kb[myMethDi_4four_sub1_10kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four_sub1_10kb[myMethDi_4four_sub1_10kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four_sub1_10kb[myMethDi_4four_sub1_10kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_4four_sub1_10kb <- -log10(myQvalue_4four_sub1_10kb)
myColor_4four_sub1_10kb[myColor_4four_sub1_10kb >  2]  = "yes"
myColor_4four_sub1_10kb[myColor_4four_sub1_10kb <= 2]  = "no"
DataFrame_Local5_4four_sub1_10kb <- data.frame(myx1 = -log10(myQvalue_4four_sub1_10kb), myy1 = myMethDi_4four_sub1_10kb,   mycolor1 = myColor_4four_sub1_10kb )

ggplot(DataFrame_Local5_4four_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_sub1_10kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_4four_sub1_10kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub1_10kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_4four_sub1_10kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_4four_sub1_10kb  = getMethylDiff(myDiff_4four_sub1_10kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four_sub1_10kb = getMethylDiff(myDiff_4four_sub1_10kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_4four_sub1_10kb       = getMethylDiff(myDiff_4four_sub1_10kb, difference=10, qvalue=0.01)
myDiffTemp_4four_sub1_10kb      = getMethylDiff(myDiff_4four_sub1_10kb, difference=0,  qvalue=0.05)

write.table(myDiff_4four_sub1_10kb , file = paste(myOutDir_4four_sub1_10kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_sub1_10kb , file = paste(myOutDir_4four_sub1_10kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_sub1_10kb , file = paste(myOutDir_4four_sub1_10kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_sub1_10kb , file = paste(myOutDir_4four_sub1_10kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_4four_sub1_10kb , file = paste(myOutDir_4four_sub1_10kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four_sub1_10kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_sub1_10kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four_sub1_10kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_sub1_10kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_4four_sub1_10kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub1_10kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub1_10kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub1_10kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffCpGann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub1_10kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffrepeatann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub1_10kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub1_10kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub1_10kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub1_10kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub1_10kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub1_10kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub1_10kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub1_10kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub1_10kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_sub1_10kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four_sub1_10kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub1_10kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub1_10kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub1_10kb,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffCpGann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub1_10kb,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffrepeatann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub1_10kb, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub1_10kb, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub1_10kb, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub1_10kb, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub1_10kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub1_10kb, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub1_10kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub1_10kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub1_10kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub1_10kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub1_10kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()







### Tiling windows analysis
##################
myOutDir_4four_sub2_5kb = paste(myOutDir_4four, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_4four_sub2_5kb) ) { dir.create(myOutDir_4four_sub2_5kb, recursive = TRUE) }


tiles_4four_sub2_5kb = tileMethylCounts(filtered.myobj_4four, win.size=5000, step.size=5000)
head(tiles_4four_sub2_5kb[[1]], 3)

meth_4four_sub2_5kb = unite(tiles_4four_sub2_5kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_4four_sub2_5kb)
dim(meth_4four_sub2_5kb)

mat_4four_sub2_5kb = percMethylation(meth_4four_sub2_5kb)
head(mat_4four_sub2_5kb)
dim(mat_4four_sub2_5kb)

write.table(meth_4four_sub2_5kb , 
            file = paste(myOutDir_4four_sub2_5kb,   "0A_meth_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub2_5kb , 
            file = paste(myOutDir_4four_sub2_5kb,   "0B_mat_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four_sub2_5kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four_sub2_5kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub2_5kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four_sub2_5kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub2_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four_sub2_5kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub2_5kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_4four_sub2_5kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four_sub2_5kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub2_5kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub2_5kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four_sub2_5kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four_sub2_5kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four_sub2_5kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub2_5kb, screeplot=TRUE)
PCASamples(meth_4four_sub2_5kb)
dev.off()



##################
myDiff_4four_sub2_5kb = calculateDiffMeth(meth_4four_sub2_5kb, num.cores=16)
dim(myDiff_4four_sub2_5kb)
names(myDiff_4four_sub2_5kb)

myQvalue_4four_sub2_5kb = myDiff_4four_sub2_5kb$qvalue
myMethDi_4four_sub2_5kb = myDiff_4four_sub2_5kb$meth.diff

pdf(paste(myOutDir_4four_sub2_5kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_4four_sub2_5kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_4four_sub2_5kb[myQvalue_4four_sub2_5kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_4four_sub2_5kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four_sub2_5kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four_sub2_5kb[myMethDi_4four_sub2_5kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four_sub2_5kb[myMethDi_4four_sub2_5kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four_sub2_5kb[myMethDi_4four_sub2_5kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_4four_sub2_5kb <- -log10(myQvalue_4four_sub2_5kb)
myColor_4four_sub2_5kb[myColor_4four_sub2_5kb >  2]  = "yes"
myColor_4four_sub2_5kb[myColor_4four_sub2_5kb <= 2]  = "no"
DataFrame_Local5_4four_sub2_5kb <- data.frame(myx1 = -log10(myQvalue_4four_sub2_5kb), myy1 = myMethDi_4four_sub2_5kb,   mycolor1 = myColor_4four_sub2_5kb )

ggplot(DataFrame_Local5_4four_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_sub2_5kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_4four_sub2_5kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub2_5kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_4four_sub2_5kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_4four_sub2_5kb  = getMethylDiff(myDiff_4four_sub2_5kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four_sub2_5kb = getMethylDiff(myDiff_4four_sub2_5kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_4four_sub2_5kb       = getMethylDiff(myDiff_4four_sub2_5kb, difference=10, qvalue=0.01)
myDiffTemp_4four_sub2_5kb      = getMethylDiff(myDiff_4four_sub2_5kb, difference=0,  qvalue=0.05)

write.table(myDiff_4four_sub2_5kb , file = paste(myOutDir_4four_sub2_5kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_sub2_5kb , file = paste(myOutDir_4four_sub2_5kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_sub2_5kb , file = paste(myOutDir_4four_sub2_5kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_sub2_5kb , file = paste(myOutDir_4four_sub2_5kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_4four_sub2_5kb , file = paste(myOutDir_4four_sub2_5kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four_sub2_5kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_sub2_5kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four_sub2_5kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_sub2_5kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_4four_sub2_5kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub2_5kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub2_5kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub2_5kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffCpGann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub2_5kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffrepeatann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub2_5kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub2_5kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub2_5kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub2_5kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub2_5kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub2_5kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub2_5kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub2_5kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub2_5kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_sub2_5kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four_sub2_5kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub2_5kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub2_5kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub2_5kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffCpGann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub2_5kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffrepeatann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub2_5kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub2_5kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub2_5kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub2_5kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub2_5kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub2_5kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub2_5kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub2_5kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub2_5kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub2_5kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub2_5kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()










### Tiling windows analysis
##################
myOutDir_4four_sub3_1kb = paste(myOutDir_4four, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_4four_sub3_1kb) ) { dir.create(myOutDir_4four_sub3_1kb, recursive = TRUE) }


tiles_4four_sub3_1kb = tileMethylCounts(filtered.myobj_4four, win.size=1000, step.size=1000)
head(tiles_4four_sub3_1kb[[1]], 3)

meth_4four_sub3_1kb = unite(tiles_4four_sub3_1kb, destrand=FALSE,  min.per.group = 5L   )
head(meth_4four_sub3_1kb)
dim(meth_4four_sub3_1kb)

mat_4four_sub3_1kb = percMethylation(meth_4four_sub3_1kb)
head(mat_4four_sub3_1kb)
dim(mat_4four_sub3_1kb)

write.table(meth_4four_sub3_1kb , 
            file = paste(myOutDir_4four_sub3_1kb,   "0A_meth_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub3_1kb , 
            file = paste(myOutDir_4four_sub3_1kb,   "0B_mat_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four_sub3_1kb, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four_sub3_1kb, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub3_1kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four_sub3_1kb, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four_sub3_1kb, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub3_1kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_4four_sub3_1kb, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four_sub3_1kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub3_1kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub3_1kb, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four_sub3_1kb, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four_sub3_1kb, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four_sub3_1kb, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub3_1kb, screeplot=TRUE)
PCASamples(meth_4four_sub3_1kb)
dev.off()



##################
myDiff_4four_sub3_1kb = calculateDiffMeth(meth_4four_sub3_1kb, num.cores=16)
dim(myDiff_4four_sub3_1kb)
names(myDiff_4four_sub3_1kb)

myQvalue_4four_sub3_1kb = myDiff_4four_sub3_1kb$qvalue
myMethDi_4four_sub3_1kb = myDiff_4four_sub3_1kb$meth.diff

pdf(paste(myOutDir_4four_sub3_1kb,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_4four_sub3_1kb, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_4four_sub3_1kb[myQvalue_4four_sub3_1kb<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_4four_sub3_1kb,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four_sub3_1kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four_sub3_1kb[myMethDi_4four_sub3_1kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four_sub3_1kb[myMethDi_4four_sub3_1kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four_sub3_1kb[myMethDi_4four_sub3_1kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_4four_sub3_1kb <- -log10(myQvalue_4four_sub3_1kb)
myColor_4four_sub3_1kb[myColor_4four_sub3_1kb >  2]  = "yes"
myColor_4four_sub3_1kb[myColor_4four_sub3_1kb <= 2]  = "no"
DataFrame_Local5_4four_sub3_1kb <- data.frame(myx1 = -log10(myQvalue_4four_sub3_1kb), myy1 = myMethDi_4four_sub3_1kb,   mycolor1 = myColor_4four_sub3_1kb )

ggplot(DataFrame_Local5_4four_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_sub3_1kb,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_4four_sub3_1kb,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_4four_sub3_1kb,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_4four_sub3_1kb  = getMethylDiff(myDiff_4four_sub3_1kb, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four_sub3_1kb = getMethylDiff(myDiff_4four_sub3_1kb, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_4four_sub3_1kb       = getMethylDiff(myDiff_4four_sub3_1kb, difference=10, qvalue=0.01)
myDiffTemp_4four_sub3_1kb      = getMethylDiff(myDiff_4four_sub3_1kb, difference=0,  qvalue=0.05)

write.table(myDiff_4four_sub3_1kb , file = paste(myOutDir_4four_sub3_1kb,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_sub3_1kb , file = paste(myOutDir_4four_sub3_1kb,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_sub3_1kb , file = paste(myOutDir_4four_sub3_1kb,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_sub3_1kb , file = paste(myOutDir_4four_sub3_1kb,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_4four_sub3_1kb , file = paste(myOutDir_4four_sub3_1kb,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four_sub3_1kb, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_sub3_1kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four_sub3_1kb, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_sub3_1kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_4four_sub3_1kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub3_1kb, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub3_1kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub3_1kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffCpGann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub3_1kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffrepeatann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub3_1kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub3_1kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub3_1kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub3_1kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub3_1kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_sub3_1kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub3_1kb, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub3_1kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub3_1kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffCpGann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub3_1kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffrepeatann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub3_1kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub3_1kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub3_1kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub3_1kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub3_1kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub3_1kb, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub3_1kb, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()













### Tiling windows analysis
##################
myOutDir_4four_sub4_500bp = paste(myOutDir_4four, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_4four_sub4_500bp) ) { dir.create(myOutDir_4four_sub4_500bp, recursive = TRUE) }


tiles_4four_sub4_500bp = tileMethylCounts(filtered.myobj_4four, win.size=500, step.size=500)
head(tiles_4four_sub4_500bp[[1]], 3)

meth_4four_sub4_500bp = unite(tiles_4four_sub4_500bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_4four_sub4_500bp)
dim(meth_4four_sub4_500bp)

mat_4four_sub4_500bp = percMethylation(meth_4four_sub4_500bp)
head(mat_4four_sub4_500bp)
dim(mat_4four_sub4_500bp)

write.table(meth_4four_sub4_500bp , 
            file = paste(myOutDir_4four_sub4_500bp,   "0A_meth_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub4_500bp , 
            file = paste(myOutDir_4four_sub4_500bp,   "0B_mat_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four_sub4_500bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four_sub4_500bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub4_500bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four_sub4_500bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four_sub4_500bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub4_500bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_4four_sub4_500bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four_sub4_500bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_500bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_500bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four_sub4_500bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four_sub4_500bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four_sub4_500bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub4_500bp, screeplot=TRUE)
PCASamples(meth_4four_sub4_500bp)
dev.off()



##################
myDiff_4four_sub4_500bp = calculateDiffMeth(meth_4four_sub4_500bp, num.cores=16)
dim(myDiff_4four_sub4_500bp)
names(myDiff_4four_sub4_500bp)

myQvalue_4four_sub4_500bp = myDiff_4four_sub4_500bp$qvalue
myMethDi_4four_sub4_500bp = myDiff_4four_sub4_500bp$meth.diff

pdf(paste(myOutDir_4four_sub4_500bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_4four_sub4_500bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_4four_sub4_500bp[myQvalue_4four_sub4_500bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_4four_sub4_500bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four_sub4_500bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_4four_sub4_500bp <- -log10(myQvalue_4four_sub4_500bp)
myColor_4four_sub4_500bp[myColor_4four_sub4_500bp >  2]  = "yes"
myColor_4four_sub4_500bp[myColor_4four_sub4_500bp <= 2]  = "no"
DataFrame_Local5_4four_sub4_500bp <- data.frame(myx1 = -log10(myQvalue_4four_sub4_500bp), myy1 = myMethDi_4four_sub4_500bp,   mycolor1 = myColor_4four_sub4_500bp )

ggplot(DataFrame_Local5_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_sub4_500bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_4four_sub4_500bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_4four_sub4_500bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_4four_sub4_500bp  = getMethylDiff(myDiff_4four_sub4_500bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four_sub4_500bp = getMethylDiff(myDiff_4four_sub4_500bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_4four_sub4_500bp       = getMethylDiff(myDiff_4four_sub4_500bp, difference=10, qvalue=0.01)
myDiffTemp_4four_sub4_500bp      = getMethylDiff(myDiff_4four_sub4_500bp, difference=0,  qvalue=0.05)

write.table(myDiff_4four_sub4_500bp , file = paste(myOutDir_4four_sub4_500bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_sub4_500bp , file = paste(myOutDir_4four_sub4_500bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_sub4_500bp , file = paste(myOutDir_4four_sub4_500bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_sub4_500bp , file = paste(myOutDir_4four_sub4_500bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_4four_sub4_500bp , file = paste(myOutDir_4four_sub4_500bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four_sub4_500bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_sub4_500bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four_sub4_500bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_sub4_500bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_4four_sub4_500bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub4_500bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub4_500bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub4_500bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffCpGann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub4_500bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffrepeatann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_sub4_500bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub4_500bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub4_500bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub4_500bp,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffCpGann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub4_500bp,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffrepeatann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub4_500bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub4_500bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
















### Tiling windows analysis
##################
myOutDir_4four_sub5_100bp = paste(myOutDir_4four, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_4four_sub5_100bp) ) { dir.create(myOutDir_4four_sub5_100bp, recursive = TRUE) }


tiles_4four_sub5_100bp = tileMethylCounts(filtered.myobj_4four, win.size=100, step.size=100)
head(tiles_4four_sub5_100bp[[1]], 3)

meth_4four_sub5_100bp = unite(tiles_4four_sub5_100bp, destrand=FALSE,  min.per.group = 5L   )
head(meth_4four_sub5_100bp)
dim(meth_4four_sub5_100bp)

mat_4four_sub5_100bp = percMethylation(meth_4four_sub5_100bp)
head(mat_4four_sub5_100bp)
dim(mat_4four_sub5_100bp)

write.table(meth_4four_sub5_100bp , 
            file = paste(myOutDir_4four_sub5_100bp,   "0A_meth_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub5_100bp , 
            file = paste(myOutDir_4four_sub5_100bp,   "0B_mat_4four_region.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four_sub5_100bp, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four_sub5_100bp, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub5_100bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four_sub5_100bp, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four_sub5_100bp, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub5_100bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_4four_sub5_100bp, "3A-clusterSamples.pdf", sep="/") , width=8, height=5  )
clusterSamples(meth_4four_sub5_100bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub5_100bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub5_100bp, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth_4four_sub5_100bp, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth_4four_sub5_100bp, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth_4four_sub5_100bp, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir_4four_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub5_100bp, screeplot=TRUE)
PCASamples(meth_4four_sub5_100bp)
dev.off()



##################
myDiff_4four_sub5_100bp = calculateDiffMeth(meth_4four_sub5_100bp, num.cores=16)
dim(myDiff_4four_sub5_100bp)
names(myDiff_4four_sub5_100bp)

myQvalue_4four_sub5_100bp = myDiff_4four_sub5_100bp$qvalue
myMethDi_4four_sub5_100bp = myDiff_4four_sub5_100bp$meth.diff

pdf(paste(myOutDir_4four_sub5_100bp,"4A_qvalue_distribution.pdf",  sep="/"))
hist(myQvalue_4four_sub5_100bp, nclass=100, xlim=c(0, 1), freq=FALSE)
hist(myQvalue_4four_sub5_100bp[myQvalue_4four_sub5_100bp<0.1], nclass=20, xlim=c(0, 0.1), freq=FALSE)
dev.off() 


pdf(paste(myOutDir_4four_sub5_100bp,  "4B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four_sub5_100bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four_sub5_100bp[myMethDi_4four_sub5_100bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four_sub5_100bp[myMethDi_4four_sub5_100bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four_sub5_100bp[myMethDi_4four_sub5_100bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
dev.off() 



myColor_4four_sub5_100bp <- -log10(myQvalue_4four_sub5_100bp)
myColor_4four_sub5_100bp[myColor_4four_sub5_100bp >  2]  = "yes"
myColor_4four_sub5_100bp[myColor_4four_sub5_100bp <= 2]  = "no"
DataFrame_Local5_4four_sub5_100bp <- data.frame(myx1 = -log10(myQvalue_4four_sub5_100bp), myy1 = myMethDi_4four_sub5_100bp,   mycolor1 = myColor_4four_sub5_100bp )

ggplot(DataFrame_Local5_4four_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_sub5_100bp,  "4C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 5) +ylim(-40, 40)
ggsave( filename = paste(myOutDir_4four_sub5_100bp,  "4C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )


ggplot(DataFrame_Local5_4four_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab("-log10(q-value)") + ylab("Difference (%)") + MyTheme_1(textSize1=14)  + xlim(0, 3) +ylim(-30, 30)
ggsave( filename = paste(myOutDir_4four_sub5_100bp,  "4C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
dev.off()





myDiff25p.hypo_4four_sub5_100bp  = getMethylDiff(myDiff_4four_sub5_100bp, difference=10, qvalue=0.01, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four_sub5_100bp = getMethylDiff(myDiff_4four_sub5_100bp, difference=10, qvalue=0.01, type="hyper")  ## more enrich in ART
myDiff25p_4four_sub5_100bp       = getMethylDiff(myDiff_4four_sub5_100bp, difference=10, qvalue=0.01)
myDiffTemp_4four_sub5_100bp      = getMethylDiff(myDiff_4four_sub5_100bp, difference=0,  qvalue=0.05)

write.table(myDiff_4four_sub5_100bp , file = paste(myOutDir_4four_sub5_100bp,"5A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_sub5_100bp , file = paste(myOutDir_4four_sub5_100bp,"5B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_sub5_100bp , file = paste(myOutDir_4four_sub5_100bp,"5C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_sub5_100bp , file = paste(myOutDir_4four_sub5_100bp,"5D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_4four_sub5_100bp , file = paste(myOutDir_4four_sub5_100bp,"5E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four_sub5_100bp, "6A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_sub5_100bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four_sub5_100bp, "6B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_sub5_100bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
dev.off()




########################## annotation for hypo sites.
diffGeneAnn_4four_sub5_100bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub5_100bp, "7A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "7B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub5_100bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub5_100bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "7C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffCpGann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "7D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub5_100bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "7E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffrepeatann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "7F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub5_100bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "8A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "8B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub5_100bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "8C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "8D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub5_100bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "8E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "8F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub5_100bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "8G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "8H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub5_100bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "8I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "8J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_sub5_100bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_sub5_100bp, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub5_100bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()

##
diffCpGann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub5_100bp,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffCpGann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()


##
diffrepeatann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub5_100bp,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffrepeatann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub5_100bp, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub5_100bp, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub5_100bp, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub5_100bp, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub5_100bp, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_sub5_100bp, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_sub5_100bp, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()
#####################################################################################################################







