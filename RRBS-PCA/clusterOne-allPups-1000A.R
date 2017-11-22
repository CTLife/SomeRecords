###########################################################################
## 1. Read the raw coverage files.
###########################################################################

myOutDir <- "1000A_rmXY_allPups_NC-vs-IVFfresh"


myFileLists <- list(
"100A-rmXY/67_NC-E24-C-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/67_NC-E24-D-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/68_NC-E56-C-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/68_NC-E56-D-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/69_NC-E123-C-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/69_NC-E123-D-Boy_Rep1.bismark.cov",
"100A-rmXY/73_NC-E114-C-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/74_NC-E54-D-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/75_NC-E98-D-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/27_E16C-boy-NC-merge_Rep1.bismark.cov",

"100A-rmXY/61_NC-BS2-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/61_NC-BS2-D-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/62_NC-BS20-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/62_NC-BS20-D-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/63_NC-E8-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/63_NC-E8-D-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/40_E12C-girl-NC_Rep1.bismark.cov",
"100A-rmXY/41_E13C-girl-NC_Rep1.bismark.cov",
"100A-rmXY/73_NC-E114-D-Girl_Rep1.bismark.cov",
"100A-rmXY/74_NC-E54-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/75_NC-E98-C-Girl-merge_Rep1.bismark.cov",

"100A-rmXY/70_ART-E113-C-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/70_ART-E113-D-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/71_ART-W58-C-Boy_Rep1.bismark.cov",
"100A-rmXY/71_ART-W58-D-Boy_Rep1.bismark.cov",
"100A-rmXY/72_ART-W779-C-Boy_Rep1.bismark.cov",
"100A-rmXY/72_W779D-ART-boy_Rep1.bismark.cov" , 
"100A-rmXY/76_ART-E18-D-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/77_ART-E69-C-Boy-merge_Rep1.bismark.cov",
"100A-rmXY/78_ART-E101-D-Boy-merge_Rep1.bismark.cov" , 
"100A-rmXY/79_ART-E72-D-Boy-merge_Rep1.bismark.cov" , 
"100A-rmXY/28_E23C-boy-IVF-fresh_Rep1.bismark.cov" , 
"100A-rmXY/29_W76C-boy-IVF-fresh_Rep1.bismark.cov" , 
"100A-rmXY/30_Q22-W1452C-boy-IVF-fresh_Rep1.bismark.cov"  ,

"100A-rmXY/64_ART-BS18-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/64_ART-BS18-D-Girl_Rep1.bismark.cov",
"100A-rmXY/65_ART-BS29-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/65_ART-BS29-D-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/66_ART-E29-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/66_ART-E29-D-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/42_W53C-girl-IVF-fresh_Rep1.bismark.cov",
"100A-rmXY/43_W81C-girl-IVF-fresh_Rep1.bismark.cov",
"100A-rmXY/44_W1694-IVF-fresh-C_Rep1.bismark.cov",
"100A-rmXY/76_ART-E18-C-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/77_ART-E69-D-Girl-merge_Rep1.bismark.cov",
"100A-rmXY/78_ART-E101-C-Girl-merge_Rep1.bismark.cov" , 
"100A-rmXY/79_ART-E72-C-Girl-merge_Rep1.bismark.cov"
)

mySampleID <- list( 
"NC1",  "NC2",   "NC3",  "NC4",  "NC5",  "NC6",  "NC7",  "NC8",   "NC9",  "NC10" ,
"NC11", "NC12",  "NC13", "NC14", "NC15", "NC16", "NC17", "NC18",  "NC19", "NC20" ,  "NC21" ,
"IVF-fresh1",  "IVF-fresh2",  "IVF-fresh3",  "IVF-fresh4" ,   "IVF-fresh5",   "IVF-fresh6",  
"IVF-fresh7",  "IVF-fresh8",  "IVF-fresh9",  "IVF-fresh10" ,  "IVF-fresh11",  "IVF-fresh12"  ,  "IVF-fresh13"  ,
"IVF-fresh14", "IVF-fresh15", "IVF-fresh16", "IVF-fresh17" ,  "IVF-fresh18",  "IVF-fresh19"  ,  "IVF-fresh20"  ,
"IVF-fresh21", "IVF-fresh22", "IVF-fresh23", "IVF-fresh24" ,  "IVF-fresh25",  "IVF-fresh26"  
)

myTreatment <- c( 
0, 0,  0, 0,  0, 0, 0, 0,  0, 0,  
0, 0,  0, 0,  0, 0, 0, 0,  0, 0, 0, 
1, 1,  1, 1,  1, 1, 1, 1,  1, 1,  1, 1, 1,
1, 1,  1, 1,  1, 1, 1, 1,  1, 1,  1, 1, 1
)       




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

meth_1one = unite( filtered.myobj_1one , destrand=FALSE,  min.per.group = 17L   )
head( meth_1one )
dim( meth_1one )

meth_2two = unite(filtered.myobj_2two, destrand=FALSE,  min.per.group = 17L   )
head(meth_2two)
dim(meth_2two)

meth_3three = unite(filtered.myobj_3three, destrand=FALSE,  min.per.group = 17L   )
head(meth_3three)
dim(meth_3three)

meth_4four = unite(filtered.myobj_4four, destrand=FALSE,  min.per.group = 17L   )
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
    clusterSamples(meth_1one, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one, screeplot=TRUE)
PCASamples(meth_1one)
dev.off()





### Tiling windows analysis
##################
myOutDir_1one_sub1_10kb = paste(myOutDir_1one, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_1one_sub1_10kb) ) { dir.create(myOutDir_1one_sub1_10kb, recursive = TRUE) }


tiles_1one_sub1_10kb = tileMethylCounts(filtered.myobj_1one, win.size=10000, step.size=10000)
head(tiles_1one_sub1_10kb[[1]], 3)

meth_1one_sub1_10kb = unite(tiles_1one_sub1_10kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub1_10kb, screeplot=TRUE)
PCASamples(meth_1one_sub1_10kb)
dev.off()






### Tiling windows analysis
##################
myOutDir_1one_sub2_5kb = paste(myOutDir_1one, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_1one_sub2_5kb) ) { dir.create(myOutDir_1one_sub2_5kb, recursive = TRUE) }


tiles_1one_sub2_5kb = tileMethylCounts(filtered.myobj_1one, win.size=5000, step.size=5000)
head(tiles_1one_sub2_5kb[[1]], 3)

meth_1one_sub2_5kb = unite(tiles_1one_sub2_5kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub2_5kb, screeplot=TRUE)
PCASamples(meth_1one_sub2_5kb)
dev.off()



### Tiling windows analysis
##################
myOutDir_1one_sub3_1kb = paste(myOutDir_1one, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_1one_sub3_1kb) ) { dir.create(myOutDir_1one_sub3_1kb, recursive = TRUE) }


tiles_1one_sub3_1kb = tileMethylCounts(filtered.myobj_1one, win.size=1000, step.size=1000)
head(tiles_1one_sub3_1kb[[1]], 3)

meth_1one_sub3_1kb = unite(tiles_1one_sub3_1kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub3_1kb, screeplot=TRUE)
PCASamples(meth_1one_sub3_1kb)
dev.off()











### Tiling windows analysis
##################
myOutDir_1one_sub4_500bp = paste(myOutDir_1one, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_1one_sub4_500bp) ) { dir.create(myOutDir_1one_sub4_500bp, recursive = TRUE) }


tiles_1one_sub4_500bp = tileMethylCounts(filtered.myobj_1one, win.size=500, step.size=500)
head(tiles_1one_sub4_500bp[[1]], 3)

meth_1one_sub4_500bp = unite(tiles_1one_sub4_500bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub4_500bp, screeplot=TRUE)
PCASamples(meth_1one_sub4_500bp)
dev.off()







### Tiling windows analysis
##################
myOutDir_1one_sub5_100bp = paste(myOutDir_1one, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_1one_sub5_100bp) ) { dir.create(myOutDir_1one_sub5_100bp, recursive = TRUE) }


tiles_1one_sub5_100bp = tileMethylCounts(filtered.myobj_1one, win.size=100, step.size=100)
head(tiles_1one_sub5_100bp[[1]], 3)

meth_1one_sub5_100bp = unite(tiles_1one_sub5_100bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_1one, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_1one, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_1one, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_1one, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_1one_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_1one_sub5_100bp, screeplot=TRUE)
PCASamples(meth_1one_sub5_100bp)
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
    clusterSamples(meth_2two, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two, screeplot=TRUE)
PCASamples(meth_2two)
dev.off()





### Tiling windows analysis
##################
myOutDir_2two_sub1_10kb = paste(myOutDir_2two, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_2two_sub1_10kb) ) { dir.create(myOutDir_2two_sub1_10kb, recursive = TRUE) }


tiles_2two_sub1_10kb = tileMethylCounts(filtered.myobj_2two, win.size=10000, step.size=10000)
head(tiles_2two_sub1_10kb[[1]], 3)

meth_2two_sub1_10kb = unite(tiles_2two_sub1_10kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub1_10kb, screeplot=TRUE)
PCASamples(meth_2two_sub1_10kb)
dev.off()






### Tiling windows analysis
##################
myOutDir_2two_sub2_5kb = paste(myOutDir_2two, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_2two_sub2_5kb) ) { dir.create(myOutDir_2two_sub2_5kb, recursive = TRUE) }


tiles_2two_sub2_5kb = tileMethylCounts(filtered.myobj_2two, win.size=5000, step.size=5000)
head(tiles_2two_sub2_5kb[[1]], 3)

meth_2two_sub2_5kb = unite(tiles_2two_sub2_5kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub2_5kb, screeplot=TRUE)
PCASamples(meth_2two_sub2_5kb)
dev.off()



### Tiling windows analysis
##################
myOutDir_2two_sub3_1kb = paste(myOutDir_2two, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_2two_sub3_1kb) ) { dir.create(myOutDir_2two_sub3_1kb, recursive = TRUE) }


tiles_2two_sub3_1kb = tileMethylCounts(filtered.myobj_2two, win.size=1000, step.size=1000)
head(tiles_2two_sub3_1kb[[1]], 3)

meth_2two_sub3_1kb = unite(tiles_2two_sub3_1kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub3_1kb, screeplot=TRUE)
PCASamples(meth_2two_sub3_1kb)
dev.off()











### Tiling windows analysis
##################
myOutDir_2two_sub4_500bp = paste(myOutDir_2two, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_2two_sub4_500bp) ) { dir.create(myOutDir_2two_sub4_500bp, recursive = TRUE) }


tiles_2two_sub4_500bp = tileMethylCounts(filtered.myobj_2two, win.size=500, step.size=500)
head(tiles_2two_sub4_500bp[[1]], 3)

meth_2two_sub4_500bp = unite(tiles_2two_sub4_500bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub4_500bp, screeplot=TRUE)
PCASamples(meth_2two_sub4_500bp)
dev.off()







### Tiling windows analysis
##################
myOutDir_2two_sub5_100bp = paste(myOutDir_2two, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_2two_sub5_100bp) ) { dir.create(myOutDir_2two_sub5_100bp, recursive = TRUE) }


tiles_2two_sub5_100bp = tileMethylCounts(filtered.myobj_2two, win.size=100, step.size=100)
head(tiles_2two_sub5_100bp[[1]], 3)

meth_2two_sub5_100bp = unite(tiles_2two_sub5_100bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_2two, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_2two, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_2two, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_2two, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_2two_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub5_100bp, screeplot=TRUE)
PCASamples(meth_2two_sub5_100bp)
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
    clusterSamples(meth_3three, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three, screeplot=TRUE)
PCASamples(meth_3three)
dev.off()





### Tiling windows analysis
##################
myOutDir_3three_sub1_10kb = paste(myOutDir_3three, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_3three_sub1_10kb) ) { dir.create(myOutDir_3three_sub1_10kb, recursive = TRUE) }


tiles_3three_sub1_10kb = tileMethylCounts(filtered.myobj_3three, win.size=10000, step.size=10000)
head(tiles_3three_sub1_10kb[[1]], 3)

meth_3three_sub1_10kb = unite(tiles_3three_sub1_10kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub1_10kb, screeplot=TRUE)
PCASamples(meth_3three_sub1_10kb)
dev.off()






### Tiling windows analysis
##################
myOutDir_3three_sub2_5kb = paste(myOutDir_3three, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_3three_sub2_5kb) ) { dir.create(myOutDir_3three_sub2_5kb, recursive = TRUE) }


tiles_3three_sub2_5kb = tileMethylCounts(filtered.myobj_3three, win.size=5000, step.size=5000)
head(tiles_3three_sub2_5kb[[1]], 3)

meth_3three_sub2_5kb = unite(tiles_3three_sub2_5kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub2_5kb, screeplot=TRUE)
PCASamples(meth_3three_sub2_5kb)
dev.off()



### Tiling windows analysis
##################
myOutDir_3three_sub3_1kb = paste(myOutDir_3three, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_3three_sub3_1kb) ) { dir.create(myOutDir_3three_sub3_1kb, recursive = TRUE) }


tiles_3three_sub3_1kb = tileMethylCounts(filtered.myobj_3three, win.size=1000, step.size=1000)
head(tiles_3three_sub3_1kb[[1]], 3)

meth_3three_sub3_1kb = unite(tiles_3three_sub3_1kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub3_1kb, screeplot=TRUE)
PCASamples(meth_3three_sub3_1kb)
dev.off()











### Tiling windows analysis
##################
myOutDir_3three_sub4_500bp = paste(myOutDir_3three, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_3three_sub4_500bp) ) { dir.create(myOutDir_3three_sub4_500bp, recursive = TRUE) }


tiles_3three_sub4_500bp = tileMethylCounts(filtered.myobj_3three, win.size=500, step.size=500)
head(tiles_3three_sub4_500bp[[1]], 3)

meth_3three_sub4_500bp = unite(tiles_3three_sub4_500bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub4_500bp, screeplot=TRUE)
PCASamples(meth_3three_sub4_500bp)
dev.off()







### Tiling windows analysis
##################
myOutDir_3three_sub5_100bp = paste(myOutDir_3three, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_3three_sub5_100bp) ) { dir.create(myOutDir_3three_sub5_100bp, recursive = TRUE) }


tiles_3three_sub5_100bp = tileMethylCounts(filtered.myobj_3three, win.size=100, step.size=100)
head(tiles_3three_sub5_100bp[[1]], 3)

meth_3three_sub5_100bp = unite(tiles_3three_sub5_100bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_3three, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_3three, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_3three, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_3three, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_3three_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub5_100bp, screeplot=TRUE)
PCASamples(meth_3three_sub5_100bp)
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
    clusterSamples(meth_4four, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four, screeplot=TRUE)
PCASamples(meth_4four)
dev.off()





### Tiling windows analysis
##################
myOutDir_4four_sub1_10kb = paste(myOutDir_4four, "/1-DMR-10kb",  sep="");
if( ! file.exists(myOutDir_4four_sub1_10kb) ) { dir.create(myOutDir_4four_sub1_10kb, recursive = TRUE) }


tiles_4four_sub1_10kb = tileMethylCounts(filtered.myobj_4four, win.size=10000, step.size=10000)
head(tiles_4four_sub1_10kb[[1]], 3)

meth_4four_sub1_10kb = unite(tiles_4four_sub1_10kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four_sub1_10kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub1_10kb, screeplot=TRUE)
PCASamples(meth_4four_sub1_10kb)
dev.off()






### Tiling windows analysis
##################
myOutDir_4four_sub2_5kb = paste(myOutDir_4four, "/2-DMR-5kb",  sep="");
if( ! file.exists(myOutDir_4four_sub2_5kb) ) { dir.create(myOutDir_4four_sub2_5kb, recursive = TRUE) }


tiles_4four_sub2_5kb = tileMethylCounts(filtered.myobj_4four, win.size=5000, step.size=5000)
head(tiles_4four_sub2_5kb[[1]], 3)

meth_4four_sub2_5kb = unite(tiles_4four_sub2_5kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four_sub2_5kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub2_5kb, screeplot=TRUE)
PCASamples(meth_4four_sub2_5kb)
dev.off()



### Tiling windows analysis
##################
myOutDir_4four_sub3_1kb = paste(myOutDir_4four, "/3-DMR-1kb",  sep="");
if( ! file.exists(myOutDir_4four_sub3_1kb) ) { dir.create(myOutDir_4four_sub3_1kb, recursive = TRUE) }


tiles_4four_sub3_1kb = tileMethylCounts(filtered.myobj_4four, win.size=1000, step.size=1000)
head(tiles_4four_sub3_1kb[[1]], 3)

meth_4four_sub3_1kb = unite(tiles_4four_sub3_1kb, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four_sub3_1kb, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub3_1kb, screeplot=TRUE)
PCASamples(meth_4four_sub3_1kb)
dev.off()











### Tiling windows analysis
##################
myOutDir_4four_sub4_500bp = paste(myOutDir_4four, "/4-DMR-500bp",  sep="");
if( ! file.exists(myOutDir_4four_sub4_500bp) ) { dir.create(myOutDir_4four_sub4_500bp, recursive = TRUE) }


tiles_4four_sub4_500bp = tileMethylCounts(filtered.myobj_4four, win.size=500, step.size=500)
head(tiles_4four_sub4_500bp[[1]], 3)

meth_4four_sub4_500bp = unite(tiles_4four_sub4_500bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four_sub4_500bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub4_500bp, screeplot=TRUE)
PCASamples(meth_4four_sub4_500bp)
dev.off()







### Tiling windows analysis
##################
myOutDir_4four_sub5_100bp = paste(myOutDir_4four, "/5-DMR-100bp",  sep="");
if( ! file.exists(myOutDir_4four_sub5_100bp) ) { dir.create(myOutDir_4four_sub5_100bp, recursive = TRUE) }


tiles_4four_sub5_100bp = tileMethylCounts(filtered.myobj_4four, win.size=100, step.size=100)
head(tiles_4four_sub5_100bp[[1]], 3)

meth_4four_sub5_100bp = unite(tiles_4four_sub5_100bp, destrand=FALSE,  min.per.group = 17L   )
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

    clusterSamples(meth_4four, dist="correlation", method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="ward",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="ward",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="single",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="single",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="complete",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="complete",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="average",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="average",     plot=TRUE)

    clusterSamples(meth_4four, dist="correlation", method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="mcquitty",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="mcquitty",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="median",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="median",     plot=TRUE)


    clusterSamples(meth_4four, dist="correlation", method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="euclidean",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="maximum",     method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="manhattan",   method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="canberra",    method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="binary",      method="centroid",     plot=TRUE)
    clusterSamples(meth_4four, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()



pdf( file=paste(myOutDir_4four_sub5_100bp, "3B-PCA.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub5_100bp, screeplot=TRUE)
PCASamples(meth_4four_sub5_100bp)
dev.off()


#####################################################################################################################


















































