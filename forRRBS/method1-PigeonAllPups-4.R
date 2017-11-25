## All reps must be overlapped. (100% overlap)
## example:  
# Rscript  method1-PigeonGirl-1.R     11A_All-Chromosomes/5_cov50reads   1011A_All-Chromosomes/5_cov50reads/method1-FemaleTwins/NC-vs-IVFfresh          



args <- commandArgs(TRUE)
print("args: ")
print(args[1])   
print(args[2])     
print("#############")

inputDir = args[1];     ## the path of input file
outDir   = args[2];     ## the path of output file
# inputDir =  "11A_All-Chromosomes/1_cov10reads"
# outDir   =  "1011A_All-Chromosomes/1_cov10reads/method1-FemaleTwins/NC-vs-IVFfresh"




###########################################################################
## 1. Read the raw coverage files.
###########################################################################
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
library(gdata)
library(ggrepel)




myFileLists <- list(
paste(inputDir, "73_NC-E114-C-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "74_NC-E54-D-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "75_NC-E98-D-Boy-merge_Rep1.bismark.cov", sep="/"),
 
paste(inputDir, "73_NC-E114-D-Girl_Rep1.bismark.cov", sep="/"),
paste(inputDir, "74_NC-E54-C-Girl-merge_Rep1.bismark.cov",       sep="/"),
paste(inputDir, "75_NC-E98-C-Girl-merge_Rep1.bismark.cov",       sep="/"),
 
paste(inputDir, "23_E45D-boy-ICSI-frozen_Rep1.bismark.cov", sep="/"),
paste(inputDir, "24_W774C-ICSI-frozen_Rep1.bismark.cov", sep="/"),
paste(inputDir, "25_E115C-boy-ICSI-frozen_Rep1.bismark.cov" ,  sep="/"),
paste(inputDir, "26_A19C-ICSI-frozen_Rep1.bismark.cov", sep="/") ,

paste(inputDir, "23_E45C-girl-ICSI-frozen_Rep1.bismark.cov", sep="/"),
paste(inputDir, "24_W774D-ICSI-frozen_Rep1.bismark.cov", sep="/"), 
paste(inputDir, "25_E115D-ICSI-frozen_Rep1.bismark.cov" , sep="/"), 
paste(inputDir, "26_A19D-ICSI-frozen_Rep1.bismark.cov" , sep="/") 

)

mySampleID <- list( 
"NC1", "NC2",  "NC3",  "NC4", "NC5",  "NC6",  
"ICSI-frozen1",  "ICSI-frozen2", "ICSI-frozen3", 
"ICSI-frozen4",  "ICSI-frozen5", "ICSI-frozen6", 
"ICSI-frozen7",  "ICSI-frozen8"  
)

myTreatment <- c( 
0, 0,  0,   0, 0,  0,  
1, 1,  1,   1, 1,  1,  1, 1 
)       

myType1 <- c(
  "73B-1", "74B-2", "75B-3",   "73G-4",   "74G-5",   "75G-6",     
  "23B-7", "24B-8", "25B-9",   "26B-10",  "23G-11",  "24G-12",   "25G-13",   "26G-14"  
)

myType2 =  c( rep(x="boy", times=3) , rep(x="girl", times=3) , rep(x="boy", times=4) ,  rep(x="girl", times=4) ) 
myType3 =  c( rep(x="NC" , times=6) , rep(x="ICSI-frozen", times=8) )

## boy=16, girl=17,  father=1, mother=2
## NC=cyan,   IVF-fresh=red, IVF-frozen=purple,  ICSI-fresh=blue, ICSI-frozen=green


myType2_shape <- c(
  "boy"=16,
  "boy"=16,
  "boy"=16,

  "girl"=17,
  "girl"=17,
  "girl"=17,

  "boy"=16,
  "boy"=16,
  "boy"=16,
  "boy"=16,

  "girl"=17,
  "girl"=17,
  "girl"=17,
  "girl"=17 
  
)


myType3_color <- c(
  "NC"="cyan",
  "NC"="cyan",
  "NC"="cyan",

  "NC"="cyan",
  "NC"="cyan",
  "NC"="cyan",
  
  "ICSI-frozen"="green",
  "ICSI-frozen"="green",
  "ICSI-frozen"="green",
  "ICSI-frozen"="green", 

  "ICSI-frozen"="green",
  "ICSI-frozen"="green",
  "ICSI-frozen"="green",
  "ICSI-frozen"="green"
)


length( myType2 )
length( myType3 )
length( myType2_shape )
length( myType3_color )




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








myOutDir <- outDir
myOutDir_sub1 = paste(myOutDir, "/1-ReadRawFiles",  sep="");
if( ! file.exists(myOutDir_sub1) ) { dir.create(myOutDir_sub1, recursive = TRUE) }

sink( file=paste(myOutDir_sub1, "1-length-variables.txt", sep="/") )
print( length(myFileLists) )
print( length(mySampleID) )
print( length(myTreatment) )
cat("############################")
print( "myTreatment:" )
print( myTreatment )
cat("############################")
print( "mySampleID:" )
print( mySampleID )
cat("############################")
print( "myFileLists:" )
print( myFileLists )
sink()


continue_on_error <- function()  {
      print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error) 



# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1, "2-theLog-of-read-inputFiles.txt", sep="/") )
myobj = methRead(myFileLists,
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


 



sink( file=paste(myOutDir_sub1, "4-dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print( "######################" )
  print(   myFileLists[[i]]  )
  print(   dim(myobj[[i]])  )
}
sink()
######################################################################################################################################################











######################################################################################################################################################
######################################################################################################################################################
myOutDir_1one = paste(myOutDir, "/2-5kbBin",  sep="");
if( ! file.exists(myOutDir_1one) ) { dir.create(myOutDir_1one, recursive = TRUE) }


#####
tiles_1one_sub1_5kb = tileMethylCounts( myobj, win.size=5000, step.size=5000) ## 5kb bin
meth_1one_sub1_5kb  = unite(tiles_1one_sub1_5kb, destrand=FALSE   )   ## 100% overlap
mat_1one_sub1_5kb   = percMethylation(meth_1one_sub1_5kb )


sink( file=paste(myOutDir_1one , "0-dimensions-5kbBin.txt", sep="/")  )
    print( tiles_1one_sub1_5kb )
    print("#########dimensions:")
    print( dim(meth_1one_sub1_5kb)  )   
    print( dim(mat_1one_sub1_5kb)   )
sink()


write.table(meth_1one_sub1_5kb , 
            file = paste(myOutDir_1one,   "1A-meth-5kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_1one_sub1_5kb , 
            file = paste(myOutDir_1one,   "1B-mat-5kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_1one, "2A-MethylationStats-5kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_1one_sub1_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_1one, "2B-MethylationStats-5kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_1one_sub1_5kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_1one, "3A-CoverageStats-5kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_1one_sub1_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_1one, "3B-CoverageStats-5kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_1one_sub1_5kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_1one, "4-clusterSamples-5kbBin.pdf", sep="/") , width=8, height=5  )

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="complete",     plot=TRUE)

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="mcquitty",     plot=TRUE)

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="median",     plot=TRUE)

    clusterSamples(meth_1one_sub1_5kb, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one_sub1_5kb, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one, "5-PCA-5kbBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub1_5kb , screeplot=TRUE)
PCASamples(meth_1one_sub1_5kb )
dev.off()

#####################


PCA_1one_sub1_5kb <- prcomp( t(mat_1one_sub1_5kb) )
names(PCA_1one_sub1_5kb)

sink( file = paste(myOutDir_1one,  "6A-PCA-5kbBin.txt",  sep="/") )
print(PCA_1one_sub1_5kb)
sink()

sink( file = paste(myOutDir_1one,  "6B-PCA-summary-5kbBin.txt",  sep="/") )
summary(PCA_1one_sub1_5kb)
sink()


sink( file = paste(myOutDir_1one,  "6C-PCA-all-5kbBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_1one_sub1_5kb$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_1one_sub1_5kb$rotation)
print("####################### myPCA2$center #########################")
print(PCA_1one_sub1_5kb$center)
print("####################### myPCA2$scale #########################")
print(PCA_1one_sub1_5kb$scale)
print("####################### myPCA2$x #########################")
print(PCA_1one_sub1_5kb$x)
sink()


pdf( file=paste(myOutDir_1one,   "6D-PCA-info-5kbBin.pdf", sep="/")  )
plot(PCA_1one_sub1_5kb, type="lines")
fviz_eig(PCA_1one_sub1_5kb)
dev.off() 




my_fviz_pca_ind1_1one_sub1_5kb <- fviz_pca_ind(PCA_1one_sub1_5kb,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_1one_sub1_5kb <- fviz_pca_ind(PCA_1one_sub1_5kb,
                                 col.ind =  as.factor( myTreatment ), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3_1one_sub1_5kb <- fviz_pca_ind(PCA_1one_sub1_5kb,
                                 col.ind =  as.factor( myTreatment ) , # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4_1one_sub1_5kb <- fviz_pca_ind(PCA_1one_sub1_5kb,
                                 col.ind =  as.factor( myTreatment ), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)


svg(file=paste(myOutDir_1one, "7A-PCA-2D-1-5kbBin.svg", sep="/") )
    print(my_fviz_pca_ind1_1one_sub1_5kb)
dev.off() 

svg(file=paste(myOutDir_1one, "7B-PCA-2D-2-5kbBin.svg", sep="/") )
    print(my_fviz_pca_ind2_1one_sub1_5kb)
dev.off() 

svg(file=paste(myOutDir_1one, "7C-PCA-2D-3-5kbBin.svg", sep="/") )
    print(my_fviz_pca_ind3_1one_sub1_5kb)
dev.off() 

svg(file=paste(myOutDir_1one, "7D-PCA-2D-4-5kbBin.svg", sep="/") )
    print(my_fviz_pca_ind4_1one_sub1_5kb)
dev.off() 



#############################


PCA_1one_sub1_5kb_matrix <- PCA_1one_sub1_5kb$x
dim( PCA_1one_sub1_5kb_matrix )

PCA_1one_sub1_5kb_Contri  <- (PCA_1one_sub1_5kb$sdev)^2
PCA_1one_sub1_5kb_Contri  <- PCA_1one_sub1_5kb_Contri/sum(PCA_1one_sub1_5kb_Contri)
PCA_1one_sub1_5kb_Contri  <- PCA_1one_sub1_5kb_Contri * 100
PCA_1one_sub1_5kb_Contri  <- round(PCA_1one_sub1_5kb_Contri, 2)

label1_1one_sub1_5kb <-   paste( "PC1 ",  "(", PCA_1one_sub1_5kb_Contri[1], "%)", sep="" )
label2_1one_sub1_5kb <-   paste( "PC2 ",  "(", PCA_1one_sub1_5kb_Contri[2], "%)", sep="" )
label3_1one_sub1_5kb <-   paste( "PC3 ",  "(", PCA_1one_sub1_5kb_Contri[3], "%)", sep="" )
label1_1one_sub1_5kb  
label2_1one_sub1_5kb  
label3_1one_sub1_5kb  


myLabel = myType1
dataframeA_1one_sub1_5kb  <- data.frame( as.data.frame(PCA_1one_sub1_5kb_matrix), myType2, myType3, myLabel   ) 
dataframeA_1one_sub1_5kb 


FigureTemp1_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_1one_sub1_5kb) +   ylab(label2_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="8A-PCA-PC1-PC2-5kbBin",  height1=4,  width1=6)


FigureTemp2_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_1one_sub1_5kb) +   ylab(label2_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="8B-PCA-PC1-PC2-alpha-5kbBin",  height1=4,  width1=6)


FigureTemp3_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_1one_sub1_5kb) +   ylab(label2_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="8C-PCA-PC1-PC2-smallDot-5kbBin",  height1=4,  width1=6)


FigureTemp4 <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_1one_sub1_5kb) +   ylab(label2_1one_sub1_5kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="8D-PCA-PC1-PC2-big-5kbBin",  height1=4,  width1=6)

 

FigureTemp5_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_1one_sub1_5kb) +   ylab(label2_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="8E-PCA-PC1-PC2-text-5kbBin",  height1=4,  width1=6)


FigureTemp6_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_1one_sub1_5kb) +   ylab(label2_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="8F-PCA-PC1-PC2-text2-5kbBin",  height1=4,  width1=6)






## plot for PC1 and PC3
FigureTemp11_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_1one_sub1_5kb) +   ylab(label3_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="9A-PCA-PC1-PC3-5kbBin",  height1=4,  width1=6)


FigureTemp12_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_1one_sub1_5kb) +   ylab(label3_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="9B-PCA-PC1-PC3-alpha-5kbBin",  height1=4,  width1=6)


FigureTemp13_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_1one_sub1_5kb) +   ylab(label3_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="9C-PCA-PC1-PC3-smallDot-5kbBin",  height1=4,  width1=6)


FigureTemp14 <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_1one_sub1_5kb) +   ylab(label3_1one_sub1_5kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="9D-PCA-PC1-PC3-big-5kbBin",  height1=4,  width1=6)



FigureTemp15_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_1one_sub1_5kb) +   ylab(label3_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="9E-PCA-PC1-PC3-text-5kbBin",  height1=4,  width1=6)


FigureTemp16_1one_sub1_5kb  <- ggplot( data = dataframeA_1one_sub1_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_1one_sub1_5kb) +   ylab(label3_1one_sub1_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_1one_sub1_5kb ,  path1=myOutDir_1one, fileName1="9F-PCA-PC1-PC3-text2-5kbBin",  height1=4,  width1=6)
##################






##################
library("scatterplot3d")

pdf( file = paste(myOutDir_1one, "10A_PCA-3d-5kbBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =myType2_shape, color=myType3_color  )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_1one, "10B_PCA-3d-5kbBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_1one, "10C_PCA-3d-label-5kbBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_1one_sub1_5kb[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_1one_sub1_5kb[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_1one_sub1_5kb[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_1one_sub1_5kb$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_1one_sub1_5kb[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()

##################################################################################################################
##################################################################################################################














######################################################################################################################################################
######################################################################################################################################################
myOutDir_2two = paste(myOutDir, "/3-1kbBin",  sep="");
if( ! file.exists(myOutDir_2two) ) { dir.create(myOutDir_2two, recursive = TRUE) }


#####
tiles_2two_sub2_1kb = tileMethylCounts( myobj, win.size=1000, step.size=1000) ## 1kb bin
meth_2two_sub2_1kb  = unite(tiles_2two_sub2_1kb, destrand=FALSE   )   ## 100% overlap
mat_2two_sub2_1kb   = percMethylation(meth_2two_sub2_1kb )


sink( file=paste(myOutDir_2two , "0-dimensions-1kbBin.txt", sep="/")  )
print( tiles_2two_sub2_1kb )
print("#########dimensions:")
print( dim(meth_2two_sub2_1kb)  )   
print( dim(mat_2two_sub2_1kb)   )
sink()


write.table(meth_2two_sub2_1kb , 
            file = paste(myOutDir_2two,   "1A-meth-1kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub2_1kb , 
            file = paste(myOutDir_2two,   "1B-mat-1kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two, "2A-MethylationStats-1kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub2_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two, "2B-MethylationStats-1kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub2_1kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two, "3A-CoverageStats-1kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub2_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two, "3B-CoverageStats-1kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub2_1kb[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_2two, "4-clusterSamples-1kbBin.pdf", sep="/") , width=8, height=5  )

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="ward",     plot=TRUE)

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="single",     plot=TRUE)

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="complete",     plot=TRUE)

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="average",     plot=TRUE)

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="mcquitty",     plot=TRUE)

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="median",     plot=TRUE)

clusterSamples(meth_2two_sub2_1kb, dist="correlation", method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="euclidean",   method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="maximum",     method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="manhattan",   method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="canberra",    method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="binary",      method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_1kb, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two, "5-PCA-1kbBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub2_1kb , screeplot=TRUE)
PCASamples(meth_2two_sub2_1kb )
dev.off()

#####################


PCA_2two_sub2_1kb <- prcomp( t(mat_2two_sub2_1kb) )
names(PCA_2two_sub2_1kb)

sink( file = paste(myOutDir_2two,  "6A-PCA-1kbBin.txt",  sep="/") )
print(PCA_2two_sub2_1kb)
sink()

sink( file = paste(myOutDir_2two,  "6B-PCA-summary-1kbBin.txt",  sep="/") )
summary(PCA_2two_sub2_1kb)
sink()


sink( file = paste(myOutDir_2two,  "6C-PCA-all-1kbBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_2two_sub2_1kb$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_2two_sub2_1kb$rotation)
print("####################### myPCA2$center #########################")
print(PCA_2two_sub2_1kb$center)
print("####################### myPCA2$scale #########################")
print(PCA_2two_sub2_1kb$scale)
print("####################### myPCA2$x #########################")
print(PCA_2two_sub2_1kb$x)
sink()


pdf( file=paste(myOutDir_2two,   "6D-PCA-info-1kbBin.pdf", sep="/")  )
plot(PCA_2two_sub2_1kb, type="lines")
fviz_eig(PCA_2two_sub2_1kb)
dev.off() 




my_fviz_pca_ind1_2two_sub2_1kb <- fviz_pca_ind(PCA_2two_sub2_1kb,
                                               col.ind = "cos2", # Color by the quality of representation
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_2two_sub2_1kb <- fviz_pca_ind(PCA_2two_sub2_1kb,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)
my_fviz_pca_ind3_2two_sub2_1kb <- fviz_pca_ind(PCA_2two_sub2_1kb,
                                               col.ind =  as.factor( myTreatment ) , # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind4_2two_sub2_1kb <- fviz_pca_ind(PCA_2two_sub2_1kb,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               #legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)


svg(file=paste(myOutDir_2two, "7A-PCA-2D-1-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_2two_sub2_1kb)
dev.off() 

svg(file=paste(myOutDir_2two, "7B-PCA-2D-2-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_2two_sub2_1kb)
dev.off() 

svg(file=paste(myOutDir_2two, "7C-PCA-2D-3-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_2two_sub2_1kb)
dev.off() 

svg(file=paste(myOutDir_2two, "7D-PCA-2D-4-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_2two_sub2_1kb)
dev.off() 



#############################


PCA_2two_sub2_1kb_matrix <- PCA_2two_sub2_1kb$x
dim( PCA_2two_sub2_1kb_matrix )

PCA_2two_sub2_1kb_Contri  <- (PCA_2two_sub2_1kb$sdev)^2
PCA_2two_sub2_1kb_Contri  <- PCA_2two_sub2_1kb_Contri/sum(PCA_2two_sub2_1kb_Contri)
PCA_2two_sub2_1kb_Contri  <- PCA_2two_sub2_1kb_Contri * 100
PCA_2two_sub2_1kb_Contri  <- round(PCA_2two_sub2_1kb_Contri, 2)

label1_2two_sub2_1kb <-   paste( "PC1 ",  "(", PCA_2two_sub2_1kb_Contri[1], "%)", sep="" )
label2_2two_sub2_1kb <-   paste( "PC2 ",  "(", PCA_2two_sub2_1kb_Contri[2], "%)", sep="" )
label3_2two_sub2_1kb <-   paste( "PC3 ",  "(", PCA_2two_sub2_1kb_Contri[3], "%)", sep="" )
label1_2two_sub2_1kb  
label2_2two_sub2_1kb  
label3_2two_sub2_1kb  


myLabel = myType1
dataframeA_2two_sub2_1kb  <- data.frame( as.data.frame(PCA_2two_sub2_1kb_matrix), myType2, myType3, myLabel   ) 
dataframeA_2two_sub2_1kb 


FigureTemp1_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_1kb) +   ylab(label2_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="8A-PCA-PC1-PC2-1kbBin",  height1=4,  width1=6)


FigureTemp2_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_1kb) +   ylab(label2_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="8B-PCA-PC1-PC2-alpha-1kbBin",  height1=4,  width1=6)


FigureTemp3_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_1kb) +   ylab(label2_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="8C-PCA-PC1-PC2-smallDot-1kbBin",  height1=4,  width1=6)


FigureTemp4 <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_1kb) +   ylab(label2_2two_sub2_1kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="8D-PCA-PC1-PC2-big-1kbBin",  height1=4,  width1=6)



FigureTemp5_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_1kb) +   ylab(label2_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="8E-PCA-PC1-PC2-text-1kbBin",  height1=4,  width1=6)


FigureTemp6_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_1kb) +   ylab(label2_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="8F-PCA-PC1-PC2-text2-1kbBin",  height1=4,  width1=6)






## plot for PC1 and PC3
FigureTemp11_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_1kb) +   ylab(label3_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="9A-PCA-PC1-PC3-1kbBin",  height1=4,  width1=6)


FigureTemp12_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_1kb) +   ylab(label3_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="9B-PCA-PC1-PC3-alpha-1kbBin",  height1=4,  width1=6)


FigureTemp13_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_1kb) +   ylab(label3_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="9C-PCA-PC1-PC3-smallDot-1kbBin",  height1=4,  width1=6)


FigureTemp14 <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_1kb) +   ylab(label3_2two_sub2_1kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="9D-PCA-PC1-PC3-big-1kbBin",  height1=4,  width1=6)



FigureTemp15_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_1kb) +   ylab(label3_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="9E-PCA-PC1-PC3-text-1kbBin",  height1=4,  width1=6)


FigureTemp16_2two_sub2_1kb  <- ggplot( data = dataframeA_2two_sub2_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_1kb) +   ylab(label3_2two_sub2_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_2two_sub2_1kb ,  path1=myOutDir_2two, fileName1="9F-PCA-PC1-PC3-text2-1kbBin",  height1=4,  width1=6)
##################






##################
library("scatterplot3d")

pdf( file = paste(myOutDir_2two, "10A_PCA-3d-1kbBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =myType2_shape, color=myType3_color, size=5 )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_2two, "10B_PCA-3d-1kbBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_2two, "10C_PCA-3d-label-1kbBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_2two_sub2_1kb[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_2two_sub2_1kb[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_2two_sub2_1kb[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_2two_sub2_1kb$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_2two_sub2_1kb[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()

##################################################################################################################
##################################################################################################################









######################################################################################################################################################
######################################################################################################################################################
myOutDir_3three = paste(myOutDir, "/4-100bpBin",  sep="");
if( ! file.exists(myOutDir_3three) ) { dir.create(myOutDir_3three, recursive = TRUE) }


#####
tiles_3three_sub3_100bp = tileMethylCounts( myobj, win.size=100, step.size=100) ## 100bp bin
meth_3three_sub3_100bp  = unite(tiles_3three_sub3_100bp, destrand=FALSE   )   ## 100% overlap
mat_3three_sub3_100bp   = percMethylation(meth_3three_sub3_100bp )


sink( file=paste(myOutDir_3three , "0-dimensions-100bpBin.txt", sep="/")  )
print( tiles_3three_sub3_100bp )
print("#########dimensions:")
print( dim(meth_3three_sub3_100bp)  )   
print( dim(mat_3three_sub3_100bp)   )
sink()


write.table(meth_3three_sub3_100bp , 
            file = paste(myOutDir_3three,   "1A-meth-100bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub3_100bp , 
            file = paste(myOutDir_3three,   "1B-mat-100bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three, "2A-MethylationStats-100bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub3_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three, "2B-MethylationStats-100bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub3_100bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three, "3A-CoverageStats-100bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub3_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three, "3B-CoverageStats-100bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub3_100bp[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_3three, "4-clusterSamples-100bpBin.pdf", sep="/") , width=8, height=5  )

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="ward",     plot=TRUE)

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="single",     plot=TRUE)

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="complete",     plot=TRUE)

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="average",     plot=TRUE)

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="mcquitty",     plot=TRUE)

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="median",     plot=TRUE)

clusterSamples(meth_3three_sub3_100bp, dist="correlation", method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="euclidean",   method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="maximum",     method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="manhattan",   method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="canberra",    method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="binary",      method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_100bp, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three, "5-PCA-100bpBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub3_100bp , screeplot=TRUE)
PCASamples(meth_3three_sub3_100bp )
dev.off()

#####################


PCA_3three_sub3_100bp <- prcomp( t(mat_3three_sub3_100bp) )
names(PCA_3three_sub3_100bp)

sink( file = paste(myOutDir_3three,  "6A-PCA-100bpBin.txt",  sep="/") )
print(PCA_3three_sub3_100bp)
sink()

sink( file = paste(myOutDir_3three,  "6B-PCA-summary-100bpBin.txt",  sep="/") )
summary(PCA_3three_sub3_100bp)
sink()


sink( file = paste(myOutDir_3three,  "6C-PCA-all-100bpBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_3three_sub3_100bp$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_3three_sub3_100bp$rotation)
print("####################### myPCA2$center #########################")
print(PCA_3three_sub3_100bp$center)
print("####################### myPCA2$scale #########################")
print(PCA_3three_sub3_100bp$scale)
print("####################### myPCA2$x #########################")
print(PCA_3three_sub3_100bp$x)
sink()


pdf( file=paste(myOutDir_3three,   "6D-PCA-info-100bpBin.pdf", sep="/")  )
plot(PCA_3three_sub3_100bp, type="lines")
fviz_eig(PCA_3three_sub3_100bp)
dev.off() 




my_fviz_pca_ind1_3three_sub3_100bp <- fviz_pca_ind(PCA_3three_sub3_100bp,
                                                   col.ind = "cos2", # Color by the quality of representation
                                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                   repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_3three_sub3_100bp <- fviz_pca_ind(PCA_3three_sub3_100bp,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   addEllipses = TRUE, # Concentration ellipses
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE
)
my_fviz_pca_ind3_3three_sub3_100bp <- fviz_pca_ind(PCA_3three_sub3_100bp,
                                                   col.ind =  as.factor( myTreatment ) , # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   addEllipses = TRUE, # Concentration ellipses
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE, 
                                                   label = "none", 
                                                   alpha.ind = 1
)
my_fviz_pca_ind4_3three_sub3_100bp <- fviz_pca_ind(PCA_3three_sub3_100bp,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   #legend.title = "Groups",
                                                   repel = TRUE, 
                                                   label = "none", 
                                                   alpha.ind = 1
)


svg(file=paste(myOutDir_3three, "7A-PCA-2D-1-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_3three_sub3_100bp)
dev.off() 

svg(file=paste(myOutDir_3three, "7B-PCA-2D-2-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_3three_sub3_100bp)
dev.off() 

svg(file=paste(myOutDir_3three, "7C-PCA-2D-3-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_3three_sub3_100bp)
dev.off() 

svg(file=paste(myOutDir_3three, "7D-PCA-2D-4-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_3three_sub3_100bp)
dev.off() 



#############################

PCA_3three_sub3_100bp_matrix <- PCA_3three_sub3_100bp$x
dim( PCA_3three_sub3_100bp_matrix )

PCA_3three_sub3_100bp_Contri  <- (PCA_3three_sub3_100bp$sdev)^2
PCA_3three_sub3_100bp_Contri  <- PCA_3three_sub3_100bp_Contri/sum(PCA_3three_sub3_100bp_Contri)
PCA_3three_sub3_100bp_Contri  <- PCA_3three_sub3_100bp_Contri * 100
PCA_3three_sub3_100bp_Contri  <- round(PCA_3three_sub3_100bp_Contri, 2)

label1_3three_sub3_100bp <-   paste( "PC1 ",  "(", PCA_3three_sub3_100bp_Contri[1], "%)", sep="" )
label2_3three_sub3_100bp <-   paste( "PC2 ",  "(", PCA_3three_sub3_100bp_Contri[2], "%)", sep="" )
label3_3three_sub3_100bp <-   paste( "PC3 ",  "(", PCA_3three_sub3_100bp_Contri[3], "%)", sep="" )
label1_3three_sub3_100bp  
label2_3three_sub3_100bp  
label3_3three_sub3_100bp  


myLabel = myType1
dataframeA_3three_sub3_100bp  <- data.frame( as.data.frame(PCA_3three_sub3_100bp_matrix), myType2, myType3, myLabel   ) 
dataframeA_3three_sub3_100bp 


FigureTemp1_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_100bp) +   ylab(label2_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="8A-PCA-PC1-PC2-100bpBin",  height1=4,  width1=6)


FigureTemp2_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_100bp) +   ylab(label2_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="8B-PCA-PC1-PC2-alpha-100bpBin",  height1=4,  width1=6)


FigureTemp3_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_100bp) +   ylab(label2_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="8C-PCA-PC1-PC2-smallDot-100bpBin",  height1=4,  width1=6)


FigureTemp4 <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_100bp) +   ylab(label2_3three_sub3_100bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="8D-PCA-PC1-PC2-big-100bpBin",  height1=4,  width1=6)



FigureTemp5_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_100bp) +   ylab(label2_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="8E-PCA-PC1-PC2-text-100bpBin",  height1=4,  width1=6)


FigureTemp6_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_100bp) +   ylab(label2_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="8F-PCA-PC1-PC2-text2-100bpBin",  height1=4,  width1=6)






## plot for PC1 and PC3
FigureTemp11_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_100bp) +   ylab(label3_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="9A-PCA-PC1-PC3-100bpBin",  height1=4,  width1=6)


FigureTemp12_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_100bp) +   ylab(label3_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="9B-PCA-PC1-PC3-alpha-100bpBin",  height1=4,  width1=6)


FigureTemp13_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_100bp) +   ylab(label3_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="9C-PCA-PC1-PC3-smallDot-100bpBin",  height1=4,  width1=6)


FigureTemp14 <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_100bp) +   ylab(label3_3three_sub3_100bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="9D-PCA-PC1-PC3-big-100bpBin",  height1=4,  width1=6)



FigureTemp15_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_100bp) +   ylab(label3_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="9E-PCA-PC1-PC3-text-100bpBin",  height1=4,  width1=6)


FigureTemp16_3three_sub3_100bp  <- ggplot( data = dataframeA_3three_sub3_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_100bp) +   ylab(label3_3three_sub3_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_3three_sub3_100bp ,  path1=myOutDir_3three, fileName1="9F-PCA-PC1-PC3-text2-100bpBin",  height1=4,  width1=6)
##################






##################
library("scatterplot3d")

pdf( file = paste(myOutDir_3three, "10A_PCA-3d-100bpBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =myType2_shape, color=myType3_color, size=5 )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_3three, "10B_PCA-3d-100bpBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_3three, "10C_PCA-3d-label-100bpBin.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_3three_sub3_100bp[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_3three_sub3_100bp[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_3three_sub3_100bp[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_3three_sub3_100bp$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_3three_sub3_100bp[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()

##################################################################################################################
##################################################################################################################




















######################################################################################################################################################
######################################################################################################################################################
myOutDir_4four = paste(myOutDir, "/5-1bp",  sep="");
if( ! file.exists(myOutDir_4four) ) { dir.create(myOutDir_4four, recursive = TRUE) }


#####
meth_4four_sub4_1bp  = unite(myobj, destrand=FALSE   )   ## 100% overlap
mat_4four_sub4_1bp   = percMethylation(meth_4four_sub4_1bp )


sink( file=paste(myOutDir_4four , "0-dimensions-1bp.txt", sep="/")  )
print(myobj)
print("#########dimensions:")
print( dim(meth_4four_sub4_1bp)  )   
print( dim(mat_4four_sub4_1bp)   )
sink()


write.table(meth_4four_sub4_1bp , 
            file = paste(myOutDir_4four,   "1A-meth-1bp.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub4_1bp , 
            file = paste(myOutDir_4four,   "1B-mat-1bp.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four, "2A-MethylationStats-1bp.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four, "2B-MethylationStats-1bp.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( myobj[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four, "3A-CoverageStats-1bp.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four, "3B-CoverageStats-1bp.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( myobj[[i]] )  )
}
sink()
sink()


pdf( file=paste(myOutDir_4four, "4-clusterSamples-1bp.pdf", sep="/") , width=8, height=5  )

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="ward",     plot=TRUE)

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="single",     plot=TRUE)

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="complete",     plot=TRUE)

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="average",     plot=TRUE)

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="mcquitty",     plot=TRUE)

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="median",     plot=TRUE)

clusterSamples(meth_4four_sub4_1bp, dist="correlation", method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="euclidean",   method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="maximum",     method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="manhattan",   method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="canberra",    method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="binary",      method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_1bp, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four, "5-PCA-1bp.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub4_1bp , screeplot=TRUE)
PCASamples(meth_4four_sub4_1bp )
dev.off()

#####################


PCA_4four_sub4_1bp <- prcomp( t(mat_4four_sub4_1bp) )
names(PCA_4four_sub4_1bp)

sink( file = paste(myOutDir_4four,  "6A-PCA-1bp.txt",  sep="/") )
print(PCA_4four_sub4_1bp)
sink()

sink( file = paste(myOutDir_4four,  "6B-PCA-summary-1bp.txt",  sep="/") )
summary(PCA_4four_sub4_1bp)
sink()


sink( file = paste(myOutDir_4four,  "6C-PCA-all-1bp.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_4four_sub4_1bp$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_4four_sub4_1bp$rotation)
print("####################### myPCA2$center #########################")
print(PCA_4four_sub4_1bp$center)
print("####################### myPCA2$scale #########################")
print(PCA_4four_sub4_1bp$scale)
print("####################### myPCA2$x #########################")
print(PCA_4four_sub4_1bp$x)
sink()


pdf( file=paste(myOutDir_4four,   "6D-PCA-info-1bp.pdf", sep="/")  )
plot(PCA_4four_sub4_1bp, type="lines")
fviz_eig(PCA_4four_sub4_1bp)
dev.off() 




my_fviz_pca_ind1_4four_sub4_1bp <- fviz_pca_ind(PCA_4four_sub4_1bp,
                                                col.ind = "cos2", # Color by the quality of representation
                                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_4four_sub4_1bp <- fviz_pca_ind(PCA_4four_sub4_1bp,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                addEllipses = TRUE, # Concentration ellipses
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE
)
my_fviz_pca_ind3_4four_sub4_1bp <- fviz_pca_ind(PCA_4four_sub4_1bp,
                                                col.ind =  as.factor( myTreatment ) , # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                addEllipses = TRUE, # Concentration ellipses
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE, 
                                                label = "none", 
                                                alpha.ind = 1
)
my_fviz_pca_ind4_4four_sub4_1bp <- fviz_pca_ind(PCA_4four_sub4_1bp,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                #legend.title = "Groups",
                                                repel = TRUE, 
                                                label = "none", 
                                                alpha.ind = 1
)


svg(file=paste(myOutDir_4four, "7A-PCA-2D-1-1bp.svg", sep="/") )
print(my_fviz_pca_ind1_4four_sub4_1bp)
dev.off() 

svg(file=paste(myOutDir_4four, "7B-PCA-2D-2-1bp.svg", sep="/") )
print(my_fviz_pca_ind2_4four_sub4_1bp)
dev.off() 

svg(file=paste(myOutDir_4four, "7C-PCA-2D-3-1bp.svg", sep="/") )
print(my_fviz_pca_ind3_4four_sub4_1bp)
dev.off() 

svg(file=paste(myOutDir_4four, "7D-PCA-2D-4-1bp.svg", sep="/") )
print(my_fviz_pca_ind4_4four_sub4_1bp)
dev.off() 



#############################
 
PCA_4four_sub4_1bp_matrix <- PCA_4four_sub4_1bp$x
dim( PCA_4four_sub4_1bp_matrix )

PCA_4four_sub4_1bp_Contri  <- (PCA_4four_sub4_1bp$sdev)^2
PCA_4four_sub4_1bp_Contri  <- PCA_4four_sub4_1bp_Contri/sum(PCA_4four_sub4_1bp_Contri)
PCA_4four_sub4_1bp_Contri  <- PCA_4four_sub4_1bp_Contri * 100
PCA_4four_sub4_1bp_Contri  <- round(PCA_4four_sub4_1bp_Contri, 2)

label1_4four_sub4_1bp <-   paste( "PC1 ",  "(", PCA_4four_sub4_1bp_Contri[1], "%)", sep="" )
label2_4four_sub4_1bp <-   paste( "PC2 ",  "(", PCA_4four_sub4_1bp_Contri[2], "%)", sep="" )
label3_4four_sub4_1bp <-   paste( "PC3 ",  "(", PCA_4four_sub4_1bp_Contri[3], "%)", sep="" )
label1_4four_sub4_1bp  
label2_4four_sub4_1bp  
label3_4four_sub4_1bp  


myLabel = myType1
dataframeA_4four_sub4_1bp  <- data.frame( as.data.frame(PCA_4four_sub4_1bp_matrix), myType2, myType3, myLabel   ) 
dataframeA_4four_sub4_1bp 


FigureTemp1_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_1bp) +   ylab(label2_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="8A-PCA-PC1-PC2-1bp",  height1=4,  width1=6)


FigureTemp2_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_1bp) +   ylab(label2_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="8B-PCA-PC1-PC2-alpha-1bp",  height1=4,  width1=6)


FigureTemp3_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_1bp) +   ylab(label2_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="8C-PCA-PC1-PC2-smallDot-1bp",  height1=4,  width1=6)


FigureTemp4 <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_1bp) +   ylab(label2_4four_sub4_1bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="8D-PCA-PC1-PC2-big-1bp",  height1=4,  width1=6)



FigureTemp5_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_1bp) +   ylab(label2_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="8E-PCA-PC1-PC2-text-1bp",  height1=4,  width1=6)


FigureTemp6_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_1bp) +   ylab(label2_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="8F-PCA-PC1-PC2-text2-1bp",  height1=4,  width1=6)






## plot for PC1 and PC3
FigureTemp11_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_1bp) +   ylab(label3_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="9A-PCA-PC1-PC3-1bp",  height1=4,  width1=6)


FigureTemp12_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_1bp) +   ylab(label3_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="9B-PCA-PC1-PC3-alpha-1bp",  height1=4,  width1=6)


FigureTemp13_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_1bp) +   ylab(label3_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="9C-PCA-PC1-PC3-smallDot-1bp",  height1=4,  width1=6)


FigureTemp14 <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_1bp) +   ylab(label3_4four_sub4_1bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="9D-PCA-PC1-PC3-big-1bp",  height1=4,  width1=6)



FigureTemp15_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_1bp) +   ylab(label3_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="9E-PCA-PC1-PC3-text-1bp",  height1=4,  width1=6)


FigureTemp16_4four_sub4_1bp  <- ggplot( data = dataframeA_4four_sub4_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_1bp) +   ylab(label3_4four_sub4_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_4four_sub4_1bp ,  path1=myOutDir_4four, fileName1="9F-PCA-PC1-PC3-text2-1bp",  height1=4,  width1=6)
##################






##################
library("scatterplot3d")

pdf( file = paste(myOutDir_4four, "10A_PCA-3d-1bp.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =myType2_shape, color=myType3_color, size=5 )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_4four, "10B_PCA-3d-1bp.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_4four, "10C_PCA-3d-label-1bp.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_4four_sub4_1bp[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_4four_sub4_1bp[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA_4four_sub4_1bp[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA_4four_sub4_1bp$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA_4four_sub4_1bp[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()

##################################################################################################################
##################################################################################################################














