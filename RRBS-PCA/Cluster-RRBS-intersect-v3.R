###########################################################################
## 1. Read the raw coverage files.
###########################################################################
myOutDir <- "2-all-CDFM-intersect"

myFileLists <- list(
"1-Coverage-CpG/1_A42C-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/1_A42D-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/1_A42F-father-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/1_A42M-Mother-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/2_W1785C-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/2_W1785D-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/2_W1785F-father-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/2_Q16-W1785M-mother-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/3_Q9-W811C-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/3_W811D-girl-IVF-frozen_Rep1.bismark.cov" ,
"1-Coverage-CpG/3_Q8-W811F-father-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/3_Q2-W811M-mother-IVF-frozen_Rep1.bismark.cov" ,  
"1-Coverage-CpG/4_W28C-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/4_W28D-girl-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/4_W28F-father-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/4_W28M-Mother-IVF-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/5_E95C-girl-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/5_E95D-girl-ICSI-fresh_Rep1.bismark.cov" , 
#"1-Coverage-CpG/5_E95D-girl-ICSI-fresh_Rep2.bismark.cov" , 
"1-Coverage-CpG/5_E95F-father-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/5_E95M-Mother-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/6_W655C-girl-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/6_W655D-girl-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/6_W655F-father-ICSI-fresh_Rep1.bismark.cov" , 
#"1-Coverage-CpG/6_W655M-Mother-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/7_W1276C-girl-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/7_W1276D-girl-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/7_W1276F-father-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/7_W1276M-mother-ICSI-fresh_Rep1.bismark.cov" , 
"1-Coverage-CpG/8_W845C-girl-ICSI-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/8_W845D-girl-ICSI-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/8_W845F-father-ICSI-frozen_Rep1.bismark.cov" , 
"1-Coverage-CpG/8_W845M-Mother-ICSI-frozen_Rep1.bismark.cov"  
)

mySampleID <- list(
"1G-1" , 
"1G-2" , 
"1F-3" , 
"1M-4" , 
"2G-5" , 
"2G-6" , 
"2F-7" , 
"2M-8" , 
"3G-9" , 
"3G-10" ,
"3F-11" , 
"3M-12" ,  
"4G-13" , 
"4G-14" , 
"4F-15" , 
"4M-16" , 
"5G-17" , 
"5G-18" , 
#"5G-19" , 
"5F-20" , 
"5M-21" , 
"6G-22" , 
"6G-23" , 
"6F-24" , 
#"6M-25" , 
"7G-26" , 
"7G-27" , 
"7F-28" , 
"7M-29" , 
"8G-30" , 
"8G-31" , 
"8F-32" , 
"8M-33" 
)

myType1 <- c(

)

myType2 <- c(
"girl" ,     
"girl" ,             
"father" ,   
"mother" ,
"girl" ,     
"girl" ,             
"father" ,   
"mother" ,
"girl" ,     
"girl" ,             
"father" ,   
"mother" ,
"girl" ,     
"girl" ,             
"father" ,   
"mother" ,
"girl" ,     
"girl" , 
#"girl" ,            
"father" ,   
"mother" ,
"girl" ,     
"girl" ,             
"father" ,   
#"mother" ,
"girl" ,     
"girl" ,             
"father" ,   
"mother" ,
"girl" ,     
"girl" ,             
"father" ,   
"mother" 
)

myType3 <- c(
"IVF-frozen" ,     
"IVF-frozen" ,           
"IVF-frozen" ,     
"IVF-frozen" ,            
"IVF-frozen" ,     
"IVF-frozen" ,             
"IVF-frozen" ,     
"IVF-frozen" ,            
"IVF-frozen" ,     
"IVF-frozen" ,            
"IVF-frozen" ,     
"IVF-frozen" ,           
"IVF-frozen" ,     
"IVF-frozen" ,              
"IVF-frozen" ,     
"IVF-frozen" ,           
"ICSI-fresh" ,     
"ICSI-fresh" ,  
#"ICSI-fresh" ,             
"ICSI-fresh" ,     
"ICSI-fresh" ,            
"ICSI-fresh" ,     
"ICSI-fresh" ,             
"ICSI-fresh" ,     
#"ICSI-fresh" ,             
"ICSI-fresh" ,     
"ICSI-fresh" ,             
"ICSI-fresh" ,     
"ICSI-fresh" ,             
"ICSI-frozen" ,     
"ICSI-frozen" ,             
"ICSI-frozen" ,     
"ICSI-frozen"  
)

## boy=16, girl=17,  father=1, mother=2
## NC=cyan,   IVF-fresh=red, IVF-frozen=purple,  ICSI-fresh=blue, ICSI-frozen=green

myType2_shape <- c(
"girl"=17 ,     
"girl"=17 ,             
"father"=1 ,   
"mother"=2 ,
"girl"=17 ,     
"girl"=17 ,             
"father"=1 ,   
"mother"=2 ,
"girl"=17 ,     
"girl"=17 ,             
"father"=1 ,   
"mother" =2,
"girl" =17,     
"girl"=17 ,             
"father"=1,   
"mother"=2 ,
"girl"=17 ,     
"girl" =17, 
#"girl"=17 ,            
"father"=1,   
"mother"=2 ,
"girl"=17 ,     
"girl" =17,             
"father"=1 ,   
#"mother"=2 ,
"girl"=17 ,     
"girl"=17 ,             
"father"=1 ,   
"mother" =2,
"girl" =17,     
"girl"=17 ,             
"father"=1 ,   
"mother"=2 
)


myType3_color <- c(
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,           
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,            
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,             
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,            
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,            
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,           
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,              
"IVF-frozen"="purple" ,     
"IVF-frozen"="purple" ,           
"ICSI-fresh"="blue" ,     
"ICSI-fresh"="blue" ,  
#"ICSI-fresh"="blue" ,             
"ICSI-fresh"="blue" ,     
"ICSI-fresh"="blue" ,            
"ICSI-fresh"="blue" ,     
"ICSI-fresh"="blue" ,             
"ICSI-fresh"="blue" ,     
#"ICSI-fresh"="blue" ,             
"ICSI-fresh"="blue" ,     
"ICSI-fresh"="blue" ,             
"ICSI-fresh"="blue" ,     
"ICSI-fresh"="blue" ,             
"ICSI-frozen"="green" ,     
"ICSI-frozen"="green" ,             
"ICSI-frozen"="green" ,     
"ICSI-frozen"="green"   
)



myOutDir_sub1 = paste(myOutDir, "/1-ReadRawFiles",  sep="");
if( ! file.exists(myOutDir_sub1) ) { dir.create(myOutDir_sub1, recursive = TRUE) }


myTreatment <- c(1 : length(myFileLists))       ## This option will determine the result of unite.
if( ! file.exists(myOutDir) ) { dir.create(myOutDir, recursive = TRUE) }

sink( file=paste(myOutDir_sub1, "1-the-Important-Variables.txt", sep="/") )
print( myTreatment )
print( length(myFileLists) )
print( length(mySampleID) )
print( length(myType2) )
print( length(myType3) )
print( length(myType2_shape) )
print( length(myType3_color) )
print( length(myTreatment) )
sink()



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
filtered.myobj_5reads = filterByCoverage(myobj,  lo.count=5,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)
sink( file=paste(myOutDir_sub1, "4B-all-files-5reads.txt", sep="/") )
    print(filtered.myobj_5reads)
sink()


filtered.myobj_10reads = filterByCoverage(myobj,  lo.count=10,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)
sink( file=paste(myOutDir_sub1, "4C-all-files-10reads.txt", sep="/") )
    print(filtered.myobj_10reads)
sink()


filtered.myobj_20reads = filterByCoverage(myobj,  lo.count=20,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)
sink( file=paste(myOutDir_sub1, "4D-all-files-20reads.txt", sep="/") )
    print(filtered.myobj_20reads)
sink()

filtered.myobj_50reads = filterByCoverage(myobj,  lo.count=50,  lo.perc=NULL,  hi.count=NULL, hi.perc=NULL)
sink( file=paste(myOutDir_sub1, "4E-all-files-50reads.txt", sep="/") )
    print(filtered.myobj_50reads)
sink()



sink( file=paste(myOutDir_sub1, "5-dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print( "######################" )
  print(   myFileLists[[i]]  )
  print(   dim(myobj[[i]])  )
  print(   dim(filtered.myobj_5reads[[i]])  )
  print(   dim(filtered.myobj_10reads[[i]])  )
  print(   dim(filtered.myobj_20reads[[i]])  )
  print(   dim(filtered.myobj_50reads[[i]])  )
}
sink()



## Merging samples
## Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. 
## This provides better coverage, but only advised when looking at CpG  
## methylation (for CpH methylation this will cause wrong results in subsequent analyses).

sink(file=paste(myOutDir_sub1, "6-log-merged-overlapSites.txt", sep="/") )

meth_1reads = unite(myobj, destrand=FALSE  )
head(meth_1reads)
dim(meth_1reads)

meth_5reads = unite(filtered.myobj_5reads, destrand=FALSE  )
head(meth_5reads)
dim(meth_5reads)

meth_10reads = unite(filtered.myobj_10reads, destrand=FALSE  )
head(meth_10reads)
dim(meth_10reads)

meth_20reads = unite(filtered.myobj_20reads, destrand=FALSE  )
head(meth_20reads)
dim(meth_20reads)

meth_50reads = unite(filtered.myobj_50reads, destrand=FALSE  )
head(meth_50reads)
dim(meth_50reads)

sink()




sink(file=paste(myOutDir_sub1, "7-dimensions-merged-overlap.txt", sep="/") )
    print( dim(meth_1reads) )  
    print( dim(meth_5reads) )
    print( dim(meth_10reads) )
    print( dim(meth_20reads) )
    print( dim(meth_50reads) )
sink()





sink(file=paste(myOutDir_sub1, "8-log-percMethylation.txt", sep="/") )

mat_1reads = percMethylation(meth_1reads)
head(mat_1reads)
dim(mat_1reads)

mat_5reads = percMethylation(meth_5reads)
head(mat_5reads)
dim(mat_5reads)

mat_10reads = percMethylation(meth_10reads)
head(mat_10reads)
dim(mat_10reads)

mat_20reads = percMethylation(meth_20reads)
head(mat_20reads)
dim(mat_20reads)

mat_50reads = percMethylation(meth_50reads)
head(mat_50reads)
dim(mat_50reads)

sink()



sink(file=paste(myOutDir_sub1, "9-dimensions-methylationMatrix.txt", sep="/") )
  print( dim(mat_1reads) )  
  print( dim(mat_5reads) )
  print( dim(mat_10reads) )
  print( dim(mat_20reads) )
  print( dim(mat_50reads) )
sink()



write.table(mat_1reads , 
            file = paste(myOutDir_sub1,"10A_mathylationLevel_1reads_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(mat_5reads , 
            file = paste(myOutDir_sub1,"10B_mathylationLevel_5reads_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(mat_10reads , 
            file = paste(myOutDir_sub1,"10C_mathylationLevel_10reads_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_20reads , 
            file = paste(myOutDir_sub1,"10D_mathylationLevel_20reads_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_50reads , 
            file = paste(myOutDir_sub1,"10E_mathylationLevel_50reads_mergeOverlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

###########################################################################








###########################################################################
## 2. All the CpG sites with 1 read at least in each sample.
###########################################################################
myOutDir_sub2 = paste(myOutDir, "/2-Cov-1reads",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

pdf( file=paste(myOutDir_sub2, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub2, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( myobj[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub2, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub2, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( myobj[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub2, "3-clusterSamples.pdf", sep="/") , width=6, height=5  )
    clusterSamples(meth_1reads, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1reads, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1reads, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_1reads, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_1reads, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_1reads, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



##################
dim(mat_1reads)
mat_1reads[1:10,1:10]
vec_1reads_values <-   array(mat_1reads) 
vec_1reads_sex    <-   rep(myType2, each = nrow(mat_1reads), times = 1)
vec_1reads_tech   <-   rep(myType3, each = nrow(mat_1reads), times = 1)
vec_1reads_fami   <-   rep( unlist(mySampleID), each = nrow(mat_1reads), times = 1)
length(vec_1reads_values)
length(vec_1reads_sex)
length(vec_1reads_tech)
length(vec_1reads_fami)
vec_1reads_values[1:10] 
vec_1reads_sex[1:10] 
vec_1reads_tech[1:10] 
vec_1reads_fami[1:10] 

DataFrame_Local2 <- data.frame( yAxis=vec_1reads_values, mySex=vec_1reads_sex, myTech=vec_1reads_tech,  myFamily=vec_1reads_fami ) 

FigureTemp2A <- ggplot(DataFrame_Local2, aes(x=myFamily ) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") +   
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,    notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="red4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
                xlab( "All samples" ) + ylab( "CpG sites mehtylation level (%)" ) + 
                ggtitle( ">= 1 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub2,  "/", "4A-violinPlot-allSamples",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=1 + length(mySampleID)/2,      dpi = 1200 )

FigureTemp2B <- ggplot(DataFrame_Local2, aes(x=mySex , y=yAxis,  fill=myTech ) ) + 
                geom_boxplot( outlier.size=0, notch=TRUE,  notchwidth = 0.1, alpha=1 ) + 
                stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="white",  geom="point", shape=19, size=3, show.legend = FALSE) + 
                xlab( "sex and age" ) + ylab( "CpG sites mehtylation level (%)" ) +  
                scale_fill_manual( values = c("NC"="cyan",   "IVF-fresh"="red", "IVF-frozen"="purple",  "ICSI-fresh"="blue", "ICSI-frozen"="green") ) +
                ggtitle( ">= 1 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub2,  "/", "4B-boxPlot-sex-tech",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=8,      dpi = 1200 )



######################################### 
myPCA2 <- prcomp( t(mat_1reads) )
names(myPCA2)

sink( file = paste(myOutDir_sub2,"5A_PCA.txt",  sep="/") )
print(myPCA2)
sink()

sink( file = paste(myOutDir_sub2,"5B_PCA-summary.txt",  sep="/") )
summary(myPCA2)
sink()

sink( file = paste(myOutDir_sub2,"5C_PCA-all.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(myPCA2$sdev)
print("####################### myPCA2$rotation #########################")
print(myPCA2$rotation)
print("####################### myPCA2$center #########################")
print(myPCA2$center)
print("####################### myPCA2$scale #########################")
print(myPCA2$scale)
print("####################### myPCA2$x #########################")
print(myPCA2$x)
sink()

pdf( file=paste(myOutDir_sub2, "5D-PCA-info.pdf", sep="/")  )
    plot(myPCA2, type="lines")
    fviz_eig(myPCA2)
dev.off() 




my_fviz_pca_ind1 <- fviz_pca_ind(myPCA2,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2 <- fviz_pca_ind(myPCA2,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3 <- fviz_pca_ind(myPCA2,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4 <- fviz_pca_ind(myPCA2,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)


svg(file=paste(myOutDir_sub2, "6A-PCA-2D-1.svg", sep="/") )
print(my_fviz_pca_ind1)
dev.off() 

svg(file=paste(myOutDir_sub2, "6B-PCA-2D-2.svg", sep="/") )
print(my_fviz_pca_ind2)
dev.off() 

svg(file=paste(myOutDir_sub2, "6C-PCA-2D-3.svg", sep="/") )
print(my_fviz_pca_ind3)
dev.off() 

svg(file=paste(myOutDir_sub2, "6D-PCA-2D-4.svg", sep="/") )
print(my_fviz_pca_ind4)
dev.off() 



myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = as.vector(unlist(mySampleID))
dataframeA  <- data.frame( as.data.frame(myPCA2_matrix), myType2, myType3, myLabel   ) 
dataframeA 

FigureTemp1 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=myOutDir_sub2, fileName1="7A-PCA-PC1-PC2",  height1=4,  width1=6)

FigureTemp2 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.7  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=myOutDir_sub2, fileName1="7B-PCA-PC1-PC2-alpha",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub2, fileName1="7C-PCA-PC1-PC2-smallDot",  height1=4,  width1=6)

FigureTemp4 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=myOutDir_sub2, fileName1="7D-PCA-PC1-PC2-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=myOutDir_sub2, fileName1="7E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6,  path1=myOutDir_sub2, fileName1="7F-PCA-PC1-PC2-text2",  height1=4,  width1=6)






####################
FigureTemp1a <- ggplot( data = dataframeA, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1a,  path1=myOutDir_sub2, fileName1="8A-PCA-PC1-PC3",  height1=4,  width1=6)

FigureTemp2a <- ggplot( data = dataframeA, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2a,  path1=myOutDir_sub2, fileName1="8B-PCA-PC1-PC3-alpha",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeA, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=1, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub2, fileName1="8C-PCA-PC1-PC3-smallDot",  height1=4,  width1=6)

FigureTemp4a <- ggplot( data = dataframeA, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4a,  path1=myOutDir_sub2, fileName1="8D-PCA-PC1-PC3-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5a <- ggplot( data = dataframeA, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5a,  path1=myOutDir_sub2, fileName1="8E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6a <- ggplot( data = dataframeA, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6a,  path1=myOutDir_sub2, fileName1="8F-PCA-PC1-PC2-text2",  height1=4,  width1=6)
                
                
                
                
                



library("scatterplot3d")

pdf( file = paste(myOutDir_sub2, "9A_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeA$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub2, "9B_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeA$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeA$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeA[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeA$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub2, "9C_PCA-3d-label.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeA[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeA$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeA$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeA[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeA$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeA[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()


#########################################################
 

 







###########################################################################
## 3. All the CpG sites with 5 read at least in each sample.
###########################################################################
myOutDir_sub3 = paste(myOutDir, "/3-Cov-5reads",  sep="");
if( ! file.exists(myOutDir_sub3) ) { dir.create(myOutDir_sub3, recursive = TRUE) }

pdf( file=paste(myOutDir_sub3, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub3, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( myobj[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub3, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub3, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( myobj[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub3, "3-clusterSamples.pdf", sep="/") , width=6, height=5  )
    clusterSamples(meth_5reads, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_5reads, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_5reads, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_5reads, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_5reads, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_5reads, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



##################
dim(mat_5reads)
mat_5reads[1:10,1:10]
vec_5reads_values <-   array(mat_5reads) 
vec_5reads_sex    <-   rep(myType2, each = nrow(mat_5reads), times = 1)
vec_5reads_tech   <-   rep(myType3, each = nrow(mat_5reads), times = 1)
vec_5reads_fami   <-   rep( unlist(mySampleID), each = nrow(mat_5reads), times = 1)
length(vec_5reads_values)
length(vec_5reads_sex)
length(vec_5reads_tech)
length(vec_5reads_fami)
vec_5reads_values[1:10] 
vec_5reads_sex[1:10] 
vec_5reads_tech[1:10] 
vec_5reads_fami[1:10] 

DataFrame_Local3 <- data.frame( yAxis=vec_5reads_values, mySex=vec_5reads_sex, myTech=vec_5reads_tech,  myFamily=vec_5reads_fami ) 

FigureTemp3A <- ggplot(DataFrame_Local3, aes(x=myFamily ) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") +   
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,    notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="red4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
                xlab( "All samples" ) + ylab( "CpG sites mehtylation level (%)" ) + 
                ggtitle( ">= 5 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub3,  "/", "4A-violinPlot-allSamples",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=1 + length(mySampleID)/2,      dpi = 1200 )

FigureTemp3B <- ggplot(DataFrame_Local3, aes(x=mySex , y=yAxis,  fill=myTech ) ) + 
                geom_boxplot( outlier.size=0, notch=TRUE,  notchwidth = 0.1, alpha=1 ) + 
                stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="white",  geom="point", shape=19, size=3, show.legend = FALSE) + 
                xlab( "sex and age" ) + ylab( "CpG sites mehtylation level (%)" ) +  
                scale_fill_manual( values = c("NC"="cyan",   "IVF-fresh"="red", "IVF-frozen"="purple",  "ICSI-fresh"="blue", "ICSI-frozen"="green") ) +
                ggtitle( ">= 5 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub3,  "/", "4B-boxPlot-sex-tech",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=8,      dpi = 1200 )



######################################### 
myPCA3 <- prcomp( t(mat_5reads) )
names(myPCA3)

sink( file = paste(myOutDir_sub3,"5A_PCA.txt",  sep="/") )
print(myPCA3)
sink()

sink( file = paste(myOutDir_sub3,"5B_PCA-summary.txt",  sep="/") )
summary(myPCA3)
sink()

sink( file = paste(myOutDir_sub3,"5C_PCA-all.txt",  sep="/") )
print("####################### myPCA3$sdev #########################")
print(myPCA3$sdev)
print("####################### myPCA3$rotation #########################")
print(myPCA3$rotation)
print("####################### myPCA3$center #########################")
print(myPCA3$center)
print("####################### myPCA3$scale #########################")
print(myPCA3$scale)
print("####################### myPCA3$x #########################")
print(myPCA3$x)
sink()

pdf( file=paste(myOutDir_sub3, "5D-PCA-info.pdf", sep="/")  )
    plot(myPCA3, type="lines")
    fviz_eig(myPCA3)
dev.off() 




my_fviz_pca_ind1 <- fviz_pca_ind(myPCA3,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2 <- fviz_pca_ind(myPCA3,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3 <- fviz_pca_ind(myPCA3,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4 <- fviz_pca_ind(myPCA3,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)


svg(file=paste(myOutDir_sub3, "6A-PCA-2D-1.svg", sep="/") )
print(my_fviz_pca_ind1)
dev.off() 

svg(file=paste(myOutDir_sub3, "6B-PCA-2D-2.svg", sep="/") )
print(my_fviz_pca_ind2)
dev.off() 

svg(file=paste(myOutDir_sub3, "6C-PCA-2D-3.svg", sep="/") )
print(my_fviz_pca_ind3)
dev.off() 

svg(file=paste(myOutDir_sub3, "6D-PCA-2D-4.svg", sep="/") )
print(my_fviz_pca_ind4)
dev.off() 



myPCA3_matrix <- myPCA3$x
dim(myPCA3_matrix)

myLabel = as.vector(unlist(mySampleID))
dataframeB  <- data.frame( as.data.frame(myPCA3_matrix), myType2, myType3, myLabel   ) 
dataframeB 

FigureTemp1 <- ggplot( data = dataframeB, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=myOutDir_sub3, fileName1="7A-PCA-PC1-PC2",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeB, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.7  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub3, fileName1="7B-PCA-PC1-PC2-alpha",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeB, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub3, fileName1="7C-PCA-PC1-PC2-smallDot",  height1=4,  width1=6)

FigureTemp4 <- ggplot( data = dataframeB, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=myOutDir_sub3, fileName1="7D-PCA-PC1-PC2-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5 <- ggplot( data = dataframeB, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=myOutDir_sub3, fileName1="7E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6 <- ggplot( data = dataframeB, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6,  path1=myOutDir_sub3, fileName1="7F-PCA-PC1-PC2-text2",  height1=4,  width1=6)






####################
FigureTemp1a <- ggplot( data = dataframeB, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1a,  path1=myOutDir_sub3, fileName1="8A-PCA-PC1-PC3",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeB, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub3, fileName1="8B-PCA-PC1-PC3-alpha",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeB, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=1, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub3, fileName1="8C-PCA-PC1-PC3-smallDot",  height1=4,  width1=6)

FigureTemp4a <- ggplot( data = dataframeB, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4a,  path1=myOutDir_sub3, fileName1="8D-PCA-PC1-PC3-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5a <- ggplot( data = dataframeB, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5a,  path1=myOutDir_sub3, fileName1="8E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6a <- ggplot( data = dataframeB, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6a,  path1=myOutDir_sub3, fileName1="8F-PCA-PC1-PC2-text2",  height1=4,  width1=6)
                
                
                
                
                



library("scatterplot3d")

pdf( file = paste(myOutDir_sub3, "9A_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeB[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeB$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeB[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeB$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeB[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeB$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub3, "9B_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeB[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeB$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeB[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeB$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeB[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeB$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub3, "9C_PCA-3d-label.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeB[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeB$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeB[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeB[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeB$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeB[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeB[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeB$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeB[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()


#########################################################
 

 







###########################################################################
## 4. All the CpG sites with 10 read at least in each sample.
###########################################################################
myOutDir_sub4 = paste(myOutDir, "/4-Cov-10reads",  sep="");
if( ! file.exists(myOutDir_sub4) ) { dir.create(myOutDir_sub4, recursive = TRUE) }

pdf( file=paste(myOutDir_sub4, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub4, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( myobj[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub4, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub4, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( myobj[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub4, "3-clusterSamples.pdf", sep="/") , width=6, height=5  )
    clusterSamples(meth_10reads, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_10reads, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_10reads, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_10reads, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_10reads, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_10reads, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



##################
dim(mat_10reads)
mat_10reads[1:10,1:10]
vec_10reads_values <-   array(mat_10reads) 
vec_10reads_sex    <-   rep(myType2, each = nrow(mat_10reads), times = 1)
vec_10reads_tech   <-   rep(myType3, each = nrow(mat_10reads), times = 1)
vec_10reads_fami   <-   rep( unlist(mySampleID), each = nrow(mat_10reads), times = 1)
length(vec_10reads_values)
length(vec_10reads_sex)
length(vec_10reads_tech)
length(vec_10reads_fami)
vec_10reads_values[1:10] 
vec_10reads_sex[1:10] 
vec_10reads_tech[1:10] 
vec_10reads_fami[1:10] 

DataFrame_Local4 <- data.frame( yAxis=vec_10reads_values, mySex=vec_10reads_sex, myTech=vec_10reads_tech,  myFamily=vec_10reads_fami ) 

FigureTemp3A <- ggplot(DataFrame_Local4, aes(x=myFamily ) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") +   
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,    notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="red4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
                xlab( "All samples" ) + ylab( "CpG sites mehtylation level (%)" ) + 
                ggtitle( ">= 10 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub4,  "/", "4A-violinPlot-allSamples",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=1 + length(mySampleID)/2,      dpi = 1200 )

FigureTemp3B <- ggplot(DataFrame_Local4, aes(x=mySex , y=yAxis,  fill=myTech ) ) + 
                geom_boxplot( outlier.size=0, notch=TRUE,  notchwidth = 0.1, alpha=1 ) + 
                stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="white",  geom="point", shape=19, size=3, show.legend = FALSE) + 
                xlab( "sex and age" ) + ylab( "CpG sites mehtylation level (%)" ) +  
                scale_fill_manual( values = c("NC"="cyan",   "IVF-fresh"="red", "IVF-frozen"="purple",  "ICSI-fresh"="blue", "ICSI-frozen"="green") ) +
                ggtitle( ">= 10 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub4,  "/", "4B-boxPlot-sex-tech",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=8,      dpi = 1200 )



######################################### 
myPCA4 <- prcomp( t(mat_10reads) )
names(myPCA4)

sink( file = paste(myOutDir_sub4,"5A_PCA.txt",  sep="/") )
print(myPCA4)
sink()

sink( file = paste(myOutDir_sub4,"5B_PCA-summary.txt",  sep="/") )
summary(myPCA4)
sink()

sink( file = paste(myOutDir_sub4,"5C_PCA-all.txt",  sep="/") )
print("####################### myPCA4$sdev #########################")
print(myPCA4$sdev)
print("####################### myPCA4$rotation #########################")
print(myPCA4$rotation)
print("####################### myPCA4$center #########################")
print(myPCA4$center)
print("####################### myPCA4$scale #########################")
print(myPCA4$scale)
print("####################### myPCA4$x #########################")
print(myPCA4$x)
sink()

pdf( file=paste(myOutDir_sub4, "5D-PCA-info.pdf", sep="/")  )
    plot(myPCA4, type="lines")
    fviz_eig(myPCA4)
dev.off() 




my_fviz_pca_ind1 <- fviz_pca_ind(myPCA4,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2 <- fviz_pca_ind(myPCA4,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3 <- fviz_pca_ind(myPCA4,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4 <- fviz_pca_ind(myPCA4,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)


svg(file=paste(myOutDir_sub4, "6A-PCA-2D-1.svg", sep="/") )
print(my_fviz_pca_ind1)
dev.off() 

svg(file=paste(myOutDir_sub4, "6B-PCA-2D-2.svg", sep="/") )
print(my_fviz_pca_ind2)
dev.off() 

svg(file=paste(myOutDir_sub4, "6C-PCA-2D-3.svg", sep="/") )
print(my_fviz_pca_ind3)
dev.off() 

svg(file=paste(myOutDir_sub4, "6D-PCA-2D-4.svg", sep="/") )
print(my_fviz_pca_ind4)
dev.off() 



myPCA4_matrix <- myPCA4$x
dim(myPCA4_matrix)

myLabel = as.vector(unlist(mySampleID))
dataframeC  <- data.frame( as.data.frame(myPCA4_matrix), myType2, myType3, myLabel   ) 
dataframeC 

FigureTemp1 <- ggplot( data = dataframeC, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=myOutDir_sub4, fileName1="7A-PCA-PC1-PC2",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeC, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.7  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub4, fileName1="7B-PCA-PC1-PC2-alpha",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeC, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub4, fileName1="7C-PCA-PC1-PC2-smallDot",  height1=4,  width1=6)

FigureTemp4 <- ggplot( data = dataframeC, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=myOutDir_sub4, fileName1="7D-PCA-PC1-PC2-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5 <- ggplot( data = dataframeC, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=myOutDir_sub4, fileName1="7E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6 <- ggplot( data = dataframeC, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6,  path1=myOutDir_sub4, fileName1="7F-PCA-PC1-PC2-text2",  height1=4,  width1=6)






####################
FigureTemp1a <- ggplot( data = dataframeC, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1a,  path1=myOutDir_sub4, fileName1="8A-PCA-PC1-PC3",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeC, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub4, fileName1="8B-PCA-PC1-PC3-alpha",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeC, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=1, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub4, fileName1="8C-PCA-PC1-PC3-smallDot",  height1=4,  width1=6)

FigureTemp4a <- ggplot( data = dataframeC, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4a,  path1=myOutDir_sub4, fileName1="8D-PCA-PC1-PC3-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5a <- ggplot( data = dataframeC, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5a,  path1=myOutDir_sub4, fileName1="8E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6a <- ggplot( data = dataframeC, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6a,  path1=myOutDir_sub4, fileName1="8F-PCA-PC1-PC2-text2",  height1=4,  width1=6)
                
                
                
                
                



library("scatterplot3d")

pdf( file = paste(myOutDir_sub4, "9A_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeC[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeC$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeC[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeC$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeC[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeC$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub4, "9B_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeC[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeC$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeC[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeC$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeC[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeC$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub4, "9C_PCA-3d-label.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeC[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeC$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeC[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeC[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeC$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeC[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeC[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeC$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeC[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()


#########################################################
 

 








###########################################################################
## 5. All the CpG sites with 20 read at least in each sample.
###########################################################################
myOutDir_sub5 = paste(myOutDir, "/5-Cov-20reads",  sep="");
if( ! file.exists(myOutDir_sub5) ) { dir.create(myOutDir_sub5, recursive = TRUE) }

pdf( file=paste(myOutDir_sub5, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub5, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( myobj[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub5, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub5, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( myobj[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub5, "3-clusterSamples.pdf", sep="/") , width=6, height=5  )
    clusterSamples(meth_20reads, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_20reads, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_20reads, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_20reads, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_20reads, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_20reads, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



##################
dim(mat_20reads)
mat_20reads[1:10,1:10]
vec_20reads_values <-   array(mat_20reads) 
vec_20reads_sex    <-   rep(myType2, each = nrow(mat_20reads), times = 1)
vec_20reads_tech   <-   rep(myType3, each = nrow(mat_20reads), times = 1)
vec_20reads_fami   <-   rep( unlist(mySampleID), each = nrow(mat_20reads), times = 1)
length(vec_20reads_values)
length(vec_20reads_sex)
length(vec_20reads_tech)
length(vec_20reads_fami)
vec_20reads_values[1:10] 
vec_20reads_sex[1:10] 
vec_20reads_tech[1:10] 
vec_20reads_fami[1:10] 

DataFrame_Local5 <- data.frame( yAxis=vec_20reads_values, mySex=vec_20reads_sex, myTech=vec_20reads_tech,  myFamily=vec_20reads_fami ) 

FigureTemp3A <- ggplot(DataFrame_Local5, aes(x=myFamily ) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") +   
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,    notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="red4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
                xlab( "All samples" ) + ylab( "CpG sites mehtylation level (%)" ) + 
                ggtitle( ">= 20 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub5,  "/", "4A-violinPlot-allSamples",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=1 + length(mySampleID)/2,      dpi = 1200 )

FigureTemp3B <- ggplot(DataFrame_Local5, aes(x=mySex , y=yAxis,  fill=myTech ) ) + 
                geom_boxplot( outlier.size=0, notch=TRUE,  notchwidth = 0.1, alpha=1 ) + 
                stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="white",  geom="point", shape=19, size=3, show.legend = FALSE) + 
                xlab( "sex and age" ) + ylab( "CpG sites mehtylation level (%)" ) +  
                scale_fill_manual( values = c("NC"="cyan",   "IVF-fresh"="red", "IVF-frozen"="purple",  "ICSI-fresh"="blue", "ICSI-frozen"="green") ) +
                ggtitle( ">= 20 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub5,  "/", "4B-boxPlot-sex-tech",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=8,      dpi = 1200 )



######################################### 
myPCA5 <- prcomp( t(mat_20reads) )
names(myPCA5)

sink( file = paste(myOutDir_sub5,"5A_PCA.txt",  sep="/") )
print(myPCA5)
sink()

sink( file = paste(myOutDir_sub5,"5B_PCA-summary.txt",  sep="/") )
summary(myPCA5)
sink()

sink( file = paste(myOutDir_sub5,"5C_PCA-all.txt",  sep="/") )
print("####################### myPCA5$sdev #########################")
print(myPCA5$sdev)
print("####################### myPCA5$rotation #########################")
print(myPCA5$rotation)
print("####################### myPCA5$center #########################")
print(myPCA5$center)
print("####################### myPCA5$scale #########################")
print(myPCA5$scale)
print("####################### myPCA5$x #########################")
print(myPCA5$x)
sink()

pdf( file=paste(myOutDir_sub5, "5D-PCA-info.pdf", sep="/")  )
    plot(myPCA5, type="lines")
    fviz_eig(myPCA5)
dev.off() 




my_fviz_pca_ind1 <- fviz_pca_ind(myPCA5,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2 <- fviz_pca_ind(myPCA5,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3 <- fviz_pca_ind(myPCA5,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4 <- fviz_pca_ind(myPCA5,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)


svg(file=paste(myOutDir_sub5, "6A-PCA-2D-1.svg", sep="/") )
print(my_fviz_pca_ind1)
dev.off() 

svg(file=paste(myOutDir_sub5, "6B-PCA-2D-2.svg", sep="/") )
print(my_fviz_pca_ind2)
dev.off() 

svg(file=paste(myOutDir_sub5, "6C-PCA-2D-3.svg", sep="/") )
print(my_fviz_pca_ind3)
dev.off() 

svg(file=paste(myOutDir_sub5, "6D-PCA-2D-4.svg", sep="/") )
print(my_fviz_pca_ind4)
dev.off() 



myPCA5_matrix <- myPCA5$x
dim(myPCA5_matrix)

myLabel = as.vector(unlist(mySampleID))
dataframeD  <- data.frame( as.data.frame(myPCA5_matrix), myType2, myType3, myLabel   ) 
dataframeD 

FigureTemp1 <- ggplot( data = dataframeD, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=myOutDir_sub5, fileName1="7A-PCA-PC1-PC2",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeD, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.7  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub5, fileName1="7B-PCA-PC1-PC2-alpha",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeD, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub5, fileName1="7C-PCA-PC1-PC2-smallDot",  height1=4,  width1=6)

FigureTemp4 <- ggplot( data = dataframeD, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=myOutDir_sub5, fileName1="7D-PCA-PC1-PC2-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5 <- ggplot( data = dataframeD, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=myOutDir_sub5, fileName1="7E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6 <- ggplot( data = dataframeD, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6,  path1=myOutDir_sub5, fileName1="7F-PCA-PC1-PC2-text2",  height1=4,  width1=6)






####################
FigureTemp1a <- ggplot( data = dataframeD, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1a,  path1=myOutDir_sub5, fileName1="8A-PCA-PC1-PC3",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeD, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub5, fileName1="8B-PCA-PC1-PC3-alpha",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeD, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=1, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub5, fileName1="8C-PCA-PC1-PC3-smallDot",  height1=4,  width1=6)

FigureTemp4a <- ggplot( data = dataframeD, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4a,  path1=myOutDir_sub5, fileName1="8D-PCA-PC1-PC3-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5a <- ggplot( data = dataframeD, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5a,  path1=myOutDir_sub5, fileName1="8E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6a <- ggplot( data = dataframeD, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6a,  path1=myOutDir_sub5, fileName1="8F-PCA-PC1-PC2-text2",  height1=4,  width1=6)
                
                
                
                
                



library("scatterplot3d")

pdf( file = paste(myOutDir_sub5, "9A_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeD[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeD$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeD[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeD$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeD[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeD$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub5, "9B_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeD[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeD$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeD[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeD$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeD[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeD$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub5, "9C_PCA-3d-label.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeD[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeD$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeD[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeD[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeD$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeD[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeD[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeD$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeD[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()


#########################################################
 

 






###########################################################################
## 6. All the CpG sites with 50 read at least in each sample.
###########################################################################
myOutDir_sub6 = paste(myOutDir, "/6-Cov-50reads",  sep="");
if( ! file.exists(myOutDir_sub6) ) { dir.create(myOutDir_sub6, recursive = TRUE) }

pdf( file=paste(myOutDir_sub6, "1A-MethylationStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub6, "1B-MethylationStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( myobj[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_sub6, "2A-CoverageStats.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_sub6, "2B-CoverageStats.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( myobj[[i]] )  )
}
sink()



pdf( file=paste(myOutDir_sub6, "3-clusterSamples.pdf", sep="/") , width=6, height=5  )
    clusterSamples(meth_50reads, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_50reads, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_50reads, dist="correlation", method="complete", plot=TRUE)
    clusterSamples(meth_50reads, dist="euclidean",   method="complete", plot=TRUE)
    clusterSamples(meth_50reads, dist="correlation", method="centroid", plot=TRUE)
    clusterSamples(meth_50reads, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



##################
dim(mat_50reads)
mat_50reads[1:10,1:10]
vec_50reads_values <-   array(mat_50reads) 
vec_50reads_sex    <-   rep(myType2, each = nrow(mat_50reads), times = 1)
vec_50reads_tech   <-   rep(myType3, each = nrow(mat_50reads), times = 1)
vec_50reads_fami   <-   rep( unlist(mySampleID), each = nrow(mat_50reads), times = 1)
length(vec_50reads_values)
length(vec_50reads_sex)
length(vec_50reads_tech)
length(vec_50reads_fami)
vec_50reads_values[1:10] 
vec_50reads_sex[1:10] 
vec_50reads_tech[1:10] 
vec_50reads_fami[1:10] 

DataFrame_Local6 <- data.frame( yAxis=vec_50reads_values, mySex=vec_50reads_sex, myTech=vec_50reads_tech,  myFamily=vec_50reads_fami ) 

FigureTemp3A <- ggplot(DataFrame_Local6, aes(x=myFamily ) ) + 
                geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") +   
                geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,    notch=FALSE,  notchwidth = 0.15, alpha=1) + 
                stat_summary( aes(y=yAxis),   fun.y=mean, colour="red4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
                xlab( "All samples" ) + ylab( "CpG sites mehtylation level (%)" ) + 
                ggtitle( ">= 20 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub6,  "/", "4A-violinPlot-allSamples",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=1 + length(mySampleID)/2,      dpi = 1200 )

FigureTemp3B <- ggplot(DataFrame_Local6, aes(x=mySex , y=yAxis,  fill=myTech ) ) + 
                geom_boxplot( outlier.size=0, notch=TRUE,  notchwidth = 0.1, alpha=1 ) + 
                stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="white",  geom="point", shape=19, size=3, show.legend = FALSE) + 
                xlab( "sex and age" ) + ylab( "CpG sites mehtylation level (%)" ) +  
                scale_fill_manual( values = c("NC"="cyan",   "IVF-fresh"="red", "IVF-frozen"="purple",  "ICSI-fresh"="blue", "ICSI-frozen"="green") ) +
                ggtitle( ">= 20 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(myOutDir_sub6,  "/", "4B-boxPlot-sex-tech",  ".svg",  sep="",  collapse=NULL),     
        height=4,    width=8,      dpi = 1200 )



######################################### 
myPCA6 <- prcomp( t(mat_50reads) )
names(myPCA6)

sink( file = paste(myOutDir_sub6,"5A_PCA.txt",  sep="/") )
print(myPCA6)
sink()

sink( file = paste(myOutDir_sub6,"5B_PCA-summary.txt",  sep="/") )
summary(myPCA6)
sink()

sink( file = paste(myOutDir_sub6,"5C_PCA-all.txt",  sep="/") )
print("####################### myPCA6$sdev #########################")
print(myPCA6$sdev)
print("####################### myPCA6$rotation #########################")
print(myPCA6$rotation)
print("####################### myPCA6$center #########################")
print(myPCA6$center)
print("####################### myPCA6$scale #########################")
print(myPCA6$scale)
print("####################### myPCA6$x #########################")
print(myPCA6$x)
sink()

pdf( file=paste(myOutDir_sub6, "5D-PCA-info.pdf", sep="/")  )
    plot(myPCA6, type="lines")
    fviz_eig(myPCA6)
dev.off() 




my_fviz_pca_ind1 <- fviz_pca_ind(myPCA6,
                                 col.ind = "cos2", # Color by the quality of representation
                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2 <- fviz_pca_ind(myPCA6,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3 <- fviz_pca_ind(myPCA6,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4 <- fviz_pca_ind(myPCA6,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)


svg(file=paste(myOutDir_sub6, "6A-PCA-2D-1.svg", sep="/") )
print(my_fviz_pca_ind1)
dev.off() 

svg(file=paste(myOutDir_sub6, "6B-PCA-2D-2.svg", sep="/") )
print(my_fviz_pca_ind2)
dev.off() 

svg(file=paste(myOutDir_sub6, "6C-PCA-2D-3.svg", sep="/") )
print(my_fviz_pca_ind3)
dev.off() 

svg(file=paste(myOutDir_sub6, "6D-PCA-2D-4.svg", sep="/") )
print(my_fviz_pca_ind4)
dev.off() 



myPCA6_matrix <- myPCA6$x
dim(myPCA6_matrix)

myLabel = as.vector(unlist(mySampleID))
dataframeE  <- data.frame( as.data.frame(myPCA6_matrix), myType2, myType3, myLabel   ) 
dataframeE 

FigureTemp1 <- ggplot( data = dataframeE, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=myOutDir_sub6, fileName1="7A-PCA-PC1-PC2",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeE, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.7  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub6, fileName1="7B-PCA-PC1-PC2-alpha",  height1=4,  width1=6)

FigureTemp3 <- ggplot( data = dataframeE, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir_sub6, fileName1="7C-PCA-PC1-PC2-smallDot",  height1=4,  width1=6)

FigureTemp4 <- ggplot( data = dataframeE, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=myOutDir_sub6, fileName1="7D-PCA-PC1-PC2-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5 <- ggplot( data = dataframeE, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=myOutDir_sub6, fileName1="7E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6 <- ggplot( data = dataframeE, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6,  path1=myOutDir_sub6, fileName1="7F-PCA-PC1-PC2-text2",  height1=4,  width1=6)






####################
FigureTemp1a <- ggplot( data = dataframeE, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1a,  path1=myOutDir_sub6, fileName1="8A-PCA-PC1-PC3",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeE, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub6, fileName1="8B-PCA-PC1-PC3-alpha",  height1=4,  width1=6)

FigureTemp3a <- ggplot( data = dataframeE, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=1, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3a,  path1=myOutDir_sub6, fileName1="8C-PCA-PC1-PC3-smallDot",  height1=4,  width1=6)

FigureTemp4a <- ggplot( data = dataframeE, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4a,  path1=myOutDir_sub6, fileName1="8D-PCA-PC1-PC3-big",  height1=4,  width1=6)

library(ggrepel)

FigureTemp5a <- ggplot( data = dataframeE, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5a,  path1=myOutDir_sub6, fileName1="8E-PCA-PC1-PC2-text",  height1=4,  width1=6)


FigureTemp6a <- ggplot( data = dataframeE, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab("PC1") +   ylab("PC3") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6a,  path1=myOutDir_sub6, fileName1="8F-PCA-PC1-PC2-text2",  height1=4,  width1=6)
                
                
                
                
                



library("scatterplot3d")

pdf( file = paste(myOutDir_sub6, "9A_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeE[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeE$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeE[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeE$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeE[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeE$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub6, "9B_PCA-3d.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeE[,1:3], pch =myType2_shape, color=myType3_color, type = "h" )
legend("top", legend = levels(dataframeE$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeE[,1:3], pch =19, color=myType3_color , type = "h")
legend("top", legend = levels(dataframeE$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

s3d <- scatterplot3d( dataframeE[,1:3], pch =myType2_shape, color="black" , type = "h")
legend("top", legend = levels(dataframeE$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)

dev.off()





pdf( file = paste(myOutDir_sub6, "9C_PCA-3d-label.pdf",  sep="/") )

s3d <- scatterplot3d( dataframeE[,1:3], pch =myType2_shape, color=myType3_color )
legend("top", legend = levels(dataframeE$myType3),
       col = unique(myType3_color), 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeE[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeE[,1:3], pch =19, color=myType3_color )
legend("top", legend = levels(dataframeE$myType3),
       col = unique(myType3_color), 
       pch = 19, 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeE[, 1:3]), labels = myLabel,  cex= 0.7 )


s3d <- scatterplot3d( dataframeE[,1:3], pch =myType2_shape, color="black" )
legend("top", legend = levels(dataframeE$myType2),
       col = "black", 
       pch = unique(myType2_shape), 
       inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d$xyz.convert(dataframeE[, 1:3]), labels = myLabel,  cex= 0.7 )


dev.off()


#########################################################
 

 







