## All reps must be overlapped. (100% overlap)
##
## example:  
## Rscript  clusterONE_FemaleTwins-1.R     11A_All-Chromosomes/5_cov50reads   1011A_All-Chromosomes/5_cov50reads/clusterONE_FemaleTwins-1.R/NC-vs-IVFfresh          



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
library(scatterplot3d)
library(car) 
library(plotly)
library(plot3D)
library(FactoMineR)
library(fpc)   



myFileLists <- list(
paste(inputDir, "67_NC-E24-C-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "67_NC-E24-D-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "68_NC-E56-C-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "68_NC-E56-D-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "69_NC-E123-C-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "69_NC-E123-D-Boy_Rep1.bismark.cov", sep="/"),
paste(inputDir, "73_NC-E114-C-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "74_NC-E54-D-Boy-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "75_NC-E98-D-Boy-merge_Rep1.bismark.cov", sep="/"),

paste(inputDir, "61_NC-BS2-C-Girl-merge_Rep1.bismark.cov",   sep="/"),
paste(inputDir, "61_NC-BS2-D-Girl-merge_Rep1.bismark.cov",   sep="/"),
paste(inputDir, "62_NC-BS20-C-Girl-merge_Rep1.bismark.cov",  sep="/"),
paste(inputDir, "62_NC-BS20-D-Girl-merge_Rep1.bismark.cov",  sep="/"),
paste(inputDir, "63_NC-E8-C-Girl_Rep3.bismark.cov",    sep="/"),
paste(inputDir, "63_NC-E8-D-Girl_Rep2.bismark.cov",    sep="/"),
paste(inputDir, "73_NC-E114-D-Girl_Rep1.bismark.cov", sep="/"),
paste(inputDir, "74_NC-E54-C-Girl-merge_Rep1.bismark.cov",       sep="/"),
paste(inputDir, "75_NC-E98-C-Girl-merge_Rep1.bismark.cov",       sep="/"),


paste(inputDir, "70_E113C-boy-ART_Rep3.bismark.cov", sep="/"),
paste(inputDir, "70_E113D-boy-ART_Rep3.bismark.cov", sep="/"),
paste(inputDir, "71_ART-W58-C-Boy_Rep1.bismark.cov", sep="/"),
paste(inputDir, "71_ART-W58-D-Boy_Rep1.bismark.cov", sep="/"),
paste(inputDir, "72_ART-W779-C-Boy_Rep1.bismark.cov", sep="/"),
paste(inputDir, "72_W779D-ART-boy_Rep1.bismark.cov",  sep="/") ,
paste(inputDir, "76_ART-E18-D-Boy_Rep2.bismark.cov", sep="/"),
paste(inputDir, "77_E69-ART-C_Rep3.bismark.cov", sep="/"),
paste(inputDir, "78_ART-E101-D-Boy_Rep4.bismark.cov" ,  sep="/"),
paste(inputDir, "79_E72-ART-D_Rep3.bismark.cov",  sep="/")  , 

paste(inputDir, "64_ART-BS18-C-Girl-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "64_ART-BS18-D-Girl_Rep1.bismark.cov",       sep="/"),
paste(inputDir, "65_ART-BS29-C-Girl-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "65_ART-BS29-D-Girl-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "66_ART-E29-C-Girl-merge_Rep1.bismark.cov",  sep="/"),
paste(inputDir, "66_ART-E29-D-Girl-merge_Rep1.bismark.cov",  sep="/")  ,
paste(inputDir, "76_ART-E18-C-Girl_Rep4.bismark.cov",    sep="/"),
paste(inputDir, "77_E69-ART-D_Rep3.bismark.cov",         sep="/"),
paste(inputDir, "78_ART-E101-C-Girl_Rep2.bismark.cov" ,  sep="/"),
paste(inputDir, "79_E72-ART-C_Rep3.bismark.cov",         sep="/") 
)




mySampleID <- list( 
"NC1",  "NC2",  "NC3",   "NC4",   "NC5",   "NC6",   "NC7",   "NC8",  "NC9",   "NC10", "NC11",
"NC12", "NC13", "NC14",  "NC15",  "NC16",  "NC17",  "NC18", 
"IVF-fresh1",  "IVF-fresh2",  "IVF-fresh3",   "IVF-fresh4",   "IVF-fresh5",   "IVF-fresh6",  
"IVF-fresh7",  "IVF-fresh8",  "IVF-fresh9",   "IVF-fresh10",  "IVF-fresh11",  "IVF-fresh12" , 
"IVF-fresh13", "IVF-fresh14", "IVF-fresh15",  "IVF-fresh16",  "IVF-fresh17",  "IVF-fresh18",  
"IVF-fresh19", "IVF-fresh20"    
)


myTreatment <- c( 
0, 0,  0, 0, 0,  0,  0, 0,  0, 0, 0,      
0, 0,  0, 0, 0,  0,  0,    
1, 1,  1, 1, 1,  1,  1, 1,  1, 1, 1,  
1, 1,  1, 1, 1,  1,  1, 1,  1 
)       


## labels of all samples
myType1 <-  as.vector( unlist(mySampleID) )


## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=red,   IVF-frozen=purple,  ICSI-fresh=blue, ICSI-frozen=green

myType2        = c( 
rep(x="boy",       times=9) ,
rep(x="girl",      times=9) ,
rep(x="boy",       times=10) ,
rep(x="girl",      times=10)  
)
myType2_shape  = c( 
rep( c("boy"=16),  times=9 ) , 
rep( c("girl"=17), times=9 ) , 
rep( c("boy"=16),  times=10 ) , 
rep( c("girl"=17), times=10 )  
) 
myType2_shape2 =  c( "boy"=16,  "girl"=17 )  



myType3        = c( rep(x="NC", times=18) , rep(x="IVF-fresh", times=20) )
myType3_color  = c( rep( c("NC"="black"), times=18) , rep( c("IVF-fresh"="red"), times=20) )
myType3_color2 = c("NC"="black" ,  "IVF-fresh"="red")

## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=red,   IVF-frozen=purple,  ICSI-fresh=blue, ICSI-frozen=green



length( myFileLists )
length( mySampleID )
length( myTreatment )
length( myType1 )
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
print( length( myFileLists ) )
print( length( mySampleID ) )
print( length( myTreatment ) )
print( length( myType1 ) )
print( length( myType2 ) )
print( length( myType3 ) )
print( length( myType2_shape ) )
print( length( myType3_color ) )
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



# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1, "2-theLog-of-read-inputFiles.txt", sep="/") )
myobj = methRead( myFileLists,
               sample.id = mySampleID,
               assembly  = "hg38",
               treatment = myTreatment,
               context   = "CpG",
               pipeline  = "bismarkCoverage",
               mincov    = 1,       ## >= n
               header    = FALSE
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



continue_on_error <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error) 

######################################################################################################################################################









######################################################################################################################################################
######################################################################################################################################################
myOutDir_2two = paste(myOutDir, "/2-500bpBin",  sep="");
if( ! file.exists(myOutDir_2two) ) { dir.create(myOutDir_2two, recursive = TRUE) }


#####
tiles_2two_sub2_500bp = tileMethylCounts( myobj, win.size=500, step.size=500) ## 500bp bin
meth_2two_sub2_500bp  = unite( tiles_2two_sub2_500bp, destrand=FALSE   )   ## 100% overlap
mat_2two_sub2_500bp   = percMethylation( meth_2two_sub2_500bp )


sink( file=paste(myOutDir_2two , "0-dimensions-500bpBin.txt", sep="/")  )
print( tiles_2two_sub2_500bp )
print("#########dimensions:")
print( dim(meth_2two_sub2_500bp)  )   
print( dim(mat_2two_sub2_500bp)   )
sink()

meth_2two_sub2_500bp[1:10,]
mat_2two_sub2_500bp[1:10,]



write.table(meth_2two_sub2_500bp , 
            file = paste(myOutDir_2two,   "1A-meth-500bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub2_500bp , 
            file = paste(myOutDir_2two,   "1B-mat-500bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_2two, "2A-MethylationStats-500bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub2_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two, "2B-MethylationStats-500bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub2_500bp[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_2two, "3A-CoverageStats-500bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub2_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two, "3B-CoverageStats-500bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub2_500bp[[i]] )  )
}
sink()
 







pdf( file=paste(myOutDir_2two, "4A-clusterSamples-default-500bpBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="ward",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="ward",     plot=TRUE)

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="single",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="single",     plot=TRUE)

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="complete",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="complete",     plot=TRUE)

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="average",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="average",     plot=TRUE)

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="mcquitty",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="mcquitty",     plot=TRUE)

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="median",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="median",     plot=TRUE)

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="centroid",     plot=TRUE)
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()








pdf( file=paste(myOutDir_2two, "4B-clusterSamples-noFilter-500bpBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()






pdf( file=paste(myOutDir_2two, "4C-clusterSamples-sd90perc-500bpBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()





pdf( file=paste(myOutDir_2two, "4D-clusterSamples-sd99perc-500bpBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()



 




pdf( file=paste(myOutDir_2two, "4E-clusterSamples-sd80perc-500bpBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_2two_sub2_500bp, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_2two_sub2_500bp, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()





pdf( file=paste(myOutDir_2two, "5-PCA-500bpBin.pdf", sep="/") , width=15, height=5  )
PCASamples(meth_2two_sub2_500bp , screeplot=TRUE)
PCASamples(meth_2two_sub2_500bp )
dev.off()

#####################





PCA_2two_sub2_500bp <- prcomp( t(mat_2two_sub2_500bp)  )
names(PCA_2two_sub2_500bp)

sink( file = paste(myOutDir_2two,  "6A-PCA-500bpBin.txt",  sep="/") )
print(PCA_2two_sub2_500bp)
sink()

sink( file = paste(myOutDir_2two,  "6B-PCA-summary-500bpBin.txt",  sep="/") )
summary(PCA_2two_sub2_500bp)
sink()


sink( file = paste(myOutDir_2two,  "6C-PCA-all-500bpBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_2two_sub2_500bp$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_2two_sub2_500bp$rotation)
print("####################### myPCA2$center #########################")
print(PCA_2two_sub2_500bp$center)
print("####################### myPCA2$scale #########################")
print(PCA_2two_sub2_500bp$scale)
print("####################### myPCA2$x #########################")
print(PCA_2two_sub2_500bp$x)
sink()


pdf( file=paste(myOutDir_2two,   "6D-PCA-info-500bpBin.pdf", sep="/")  )
plot(PCA_2two_sub2_500bp, type="lines")
fviz_eig(PCA_2two_sub2_500bp)
dev.off() 




my_fviz_pca_ind1_2two_sub2_500bp <- fviz_pca_ind(PCA_2two_sub2_500bp,
                                               col.ind = "cos2", # Color by the quality of representation
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_2two_sub2_500bp <- fviz_pca_ind(PCA_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)
my_fviz_pca_ind3_2two_sub2_500bp <- fviz_pca_ind(PCA_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ) , # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind4_2two_sub2_500bp <- fviz_pca_ind(PCA_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               #legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)


svg(file=paste(myOutDir_2two, "7A-PCA-2D-1-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "7B-PCA-2D-2-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "7C-PCA-2D-3-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "7D-PCA-2D-4-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_2two_sub2_500bp)
dev.off() 



#############################


PCA_2two_sub2_500bp_matrix <- PCA_2two_sub2_500bp$x
dim( PCA_2two_sub2_500bp_matrix )

PCA_2two_sub2_500bp_Contri  <- (PCA_2two_sub2_500bp$sdev)^2
PCA_2two_sub2_500bp_Contri  <- PCA_2two_sub2_500bp_Contri/sum(PCA_2two_sub2_500bp_Contri)
PCA_2two_sub2_500bp_Contri  <- PCA_2two_sub2_500bp_Contri * 100
PCA_2two_sub2_500bp_Contri  <- round(PCA_2two_sub2_500bp_Contri, 2)

label1_2two_sub2_500bp <-   paste( "PC1 ",  "(", PCA_2two_sub2_500bp_Contri[1], "%)", sep="" )
label2_2two_sub2_500bp <-   paste( "PC2 ",  "(", PCA_2two_sub2_500bp_Contri[2], "%)", sep="" )
label3_2two_sub2_500bp <-   paste( "PC3 ",  "(", PCA_2two_sub2_500bp_Contri[3], "%)", sep="" )
label1_2two_sub2_500bp  
label2_2two_sub2_500bp  
label3_2two_sub2_500bp  


myLabel = myType1
dataframeA_2two_sub2_500bp  <- data.frame( as.data.frame(PCA_2two_sub2_500bp_matrix), myType2, myType3, myLabel   ) 
dataframeA_2two_sub2_500bp 


FigureTemp1_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="8A-PCA-PC1-PC2-500bpBin",  height1=3,  width1=4.8)


FigureTemp2_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="8B-PCA-PC1-PC2-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp3_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="8C-PCA-PC1-PC2-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="8D-PCA-PC1-PC2-big-500bpBin",   height1=3,  width1=4.8)



FigureTemp5_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="8E-PCA-PC1-PC2-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp6_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="8F-PCA-PC1-PC2-text2-500bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="9A-PCA-PC1-PC3-500bpBin",   height1=3,  width1=4.8)


FigureTemp12_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="9B-PCA-PC1-PC3-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp13_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="9C-PCA-PC1-PC3-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="9D-PCA-PC1-PC3-big-500bpBin",  height1=3,  width1=4.8)



FigureTemp15_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="9E-PCA-PC1-PC3-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp16_2two_sub2_500bp  <- ggplot( data = dataframeA_2two_sub2_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="9F-PCA-PC1-PC3-text2-500bpBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_2two, "10A_PCA-3d-500bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_2two_sub2_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
  )
legend("top", legend = levels(dataframeA_2two_sub2_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_2two_sub2_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_2two_sub2_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_2two_sub2_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_2two_sub2_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_2two, "10B_PCA-3d-500bpBin-withLabels.pdf",  sep="/") )

s3d1_2two_sub2_500bp <- scatterplot3d( 
  dataframeA_2two_sub2_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_2two_sub2_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_2two_sub2_500bp$xyz.convert(dataframeA_2two_sub2_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_2two_sub2_500bp <- scatterplot3d( 
  dataframeA_2two_sub2_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_2two_sub2_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_2two_sub2_500bp$xyz.convert(dataframeA_2two_sub2_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_2two_sub2_500bp <- scatterplot3d( 
  dataframeA_2two_sub2_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_2two_sub2_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_2two_sub2_500bp$xyz.convert(dataframeA_2two_sub2_500bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_2two, "11A_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_2two, "11B_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_2two, "12_PCA-3d-500bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )

 

##
scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )



##
scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeA_2two_sub2_500bp[,1], y = dataframeA_2two_sub2_500bp[,2], z = dataframeA_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.

 
myres.pca_2two_sub2_500bp  <- PCA(t(mat_2two_sub2_500bp),   graph = FALSE)

sink( file = paste(myOutDir_2two, "20_PCA-info-500bpBin-byPCA.txt",  sep="/") )
print( myres.pca_2two_sub2_500bp )
print( "#################################" )
print( summary(myres.pca_2two_sub2_500bp) )
print( "#################################"  )
myeig.val_2two_sub2_500bp <- get_eigenvalue( myres.pca_2two_sub2_500bp )
myeig.val_2two_sub2_500bp
sink() 


pdf( file = paste(myOutDir_2two, "21_PCA-screePlot-500bpBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_2two_sub2_500bp, addlabels = TRUE )
fviz_screeplot(X=myres.pca_2two_sub2_500bp, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_2two_sub2_500bp, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_2two_sub2_500bp, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_2two_sub2_500bp, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_2two_sub2_500bp <- fviz_pca_ind(myres.pca_2two_sub2_500bp,
                                               col.ind = "cos2", # Color by the quality of representation
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_2two_sub2_500bp <- fviz_pca_ind(myres.pca_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)
my_fviz_pca_ind3_2two_sub2_500bp <- fviz_pca_ind(myres.pca_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ) , # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind4_2two_sub2_500bp <- fviz_pca_ind(myres.pca_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               #legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind5_2two_sub2_500bp <- fviz_pca_ind(myres.pca_2two_sub2_500bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)

svg(file=paste(myOutDir_2two, "22A-PCA-2D-1-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "22B-PCA-2D-2-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "22C-PCA-2D-3-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "22D-PCA-2D-4-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_2two_sub2_500bp)
dev.off() 

svg(file=paste(myOutDir_2two, "22E-PCA-2D-4-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind5_2two_sub2_500bp)
dev.off() 



#############################


myres.pca_2two_sub2_500bp_matrix <- myres.pca_2two_sub2_500bp$ind$coord
dim( myres.pca_2two_sub2_500bp_matrix )

myres.pca_2two_sub2_500bp_Contri  <- (myres.pca_2two_sub2_500bp$eig)[,2] 
myres.pca_2two_sub2_500bp_Contri  <- myres.pca_2two_sub2_500bp_Contri/sum(myres.pca_2two_sub2_500bp_Contri)
myres.pca_2two_sub2_500bp_Contri  <- myres.pca_2two_sub2_500bp_Contri * 100
myres.pca_2two_sub2_500bp_Contri  <- round(myres.pca_2two_sub2_500bp_Contri, 2)

label1_2two_sub2_500bp <-   paste( "PC1 ",  "(", myres.pca_2two_sub2_500bp_Contri[1], "%)", sep="" )
label2_2two_sub2_500bp <-   paste( "PC2 ",  "(", myres.pca_2two_sub2_500bp_Contri[2], "%)", sep="" )
label3_2two_sub2_500bp <-   paste( "PC3 ",  "(", myres.pca_2two_sub2_500bp_Contri[3], "%)", sep="" )
label1_2two_sub2_500bp  
label2_2two_sub2_500bp  
label3_2two_sub2_500bp  


myLabel = myType1
dataframeB_2two_sub2_500bp  <- data.frame( as.data.frame(myres.pca_2two_sub2_500bp_matrix), myType2, myType3, myLabel   ) 
dataframeB_2two_sub2_500bp 


FigureTemp1_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="23A-PCA-PC1-PC2-500bpBin",  height1=3,  width1=4.8)


FigureTemp2_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="23B-PCA-PC1-PC2-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp3_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="23C-PCA-PC1-PC2-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="23D-PCA-PC1-PC2-big-500bpBin",   height1=3,  width1=4.8)



FigureTemp5_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="23E-PCA-PC1-PC2-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp6_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label2_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="23F-PCA-PC1-PC2-text2-500bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="24A-PCA-PC1-PC3-500bpBin",   height1=3,  width1=4.8)


FigureTemp12_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="24B-PCA-PC1-PC3-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp13_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="24C-PCA-PC1-PC3-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="24D-PCA-PC1-PC3-big-500bpBin",  height1=3,  width1=4.8)



FigureTemp15_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="24E-PCA-PC1-PC3-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp16_2two_sub2_500bp  <- ggplot( data = dataframeB_2two_sub2_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_500bp) +   ylab(label3_2two_sub2_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_2two_sub2_500bp ,  path1=myOutDir_2two, fileName1="24F-PCA-PC1-PC3-text2-500bpBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_2two, "25A_PCA-3d-500bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_2two_sub2_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_2two_sub2_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_2two_sub2_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_2two, "25B_PCA-3d-500bpBin-withLabels.pdf",  sep="/") )

s3d1_2two_sub2_500bp <- scatterplot3d( 
  dataframeB_2two_sub2_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_2two_sub2_500bp$xyz.convert(dataframeB_2two_sub2_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_2two_sub2_500bp <- scatterplot3d( 
  dataframeB_2two_sub2_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_2two_sub2_500bp$xyz.convert(dataframeB_2two_sub2_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_2two_sub2_500bp <- scatterplot3d( 
  dataframeB_2two_sub2_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_2two_sub2_500bp$xyz.convert(dataframeB_2two_sub2_500bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_2two, "26A_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_2two, "26B_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_2two, "27_PCA-3d-500bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )



##
scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )



##
scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )


scatter3D(x = dataframeB_2two_sub2_500bp[,1], y = dataframeB_2two_sub2_500bp[,2], z = dataframeB_2two_sub2_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_500bp,  ylab = label2_2two_sub2_500bp,  zlab = label3_2two_sub2_500bp  )

dev.off()






#################



km.res1_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=2, iter.max = 20 )
km.res2_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=3, iter.max = 20  )
km.res3_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=4, iter.max = 20  )
km.res4_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=5, iter.max = 20 )
km.res5_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=6, iter.max = 20 )
km.res6_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=7, iter.max = 20 )
km.res7_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=8, iter.max = 20 )
km.res8_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=9, iter.max = 20 )
km.res9_2two_sub2_500bp  <- kmeans( t(mat_2two_sub2_500bp), centers=10, iter.max = 20 )
km.res10_2two_sub2_500bp <- kmeans( t(mat_2two_sub2_500bp), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_2two, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_2two_sub2_500bp$cluster
print("#################### 3 classes:")
km.res2_2two_sub2_500bp$cluster
print("#################### 4 classes:")
km.res3_2two_sub2_500bp$cluster
print("#################### 5 classes:")
km.res4_2two_sub2_500bp$cluster
print("#################### 6 classes:")
km.res5_2two_sub2_500bp$cluster
print("#################### 7 classes:")
km.res6_2two_sub2_500bp$cluster
print("#################### 8 classes:")
km.res7_2two_sub2_500bp$cluster
print("#################### 9 classes:")
km.res8_2two_sub2_500bp$cluster
print("#################### 10 classes:")
km.res9_2two_sub2_500bp$cluster
print("#################### 11 classes:")
km.res10_2two_sub2_500bp$cluster
sink()



pdf( file = paste(myOutDir_2two, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_500bp), col = km.res1_2two_sub2_500bp$cluster )

plot( t(mat_2two_sub2_500bp), col = km.res1_2two_sub2_500bp$cluster)
points(km.res1_2two_sub2_500bp$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_500bp), col = km.res1_2two_sub2_500bp$cluster)
points(km.res1_2two_sub2_500bp$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_2two_sub2_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 
 

 



pdf( file = paste(myOutDir_2two, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_500bp), col = km.res2_2two_sub2_500bp$cluster )

plot( t(mat_2two_sub2_500bp), col = km.res2_2two_sub2_500bp$cluster)
points(km.res2_2two_sub2_500bp$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_500bp), col = km.res2_2two_sub2_500bp$cluster)
points(km.res2_2two_sub2_500bp$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_2two_sub2_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_2two, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_500bp), col = km.res3_2two_sub2_500bp$cluster )

plot( t(mat_2two_sub2_500bp), col = km.res3_2two_sub2_500bp$cluster)
points(km.res3_2two_sub2_500bp$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_500bp), col = km.res3_2two_sub2_500bp$cluster)
points(km.res3_2two_sub2_500bp$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_2two_sub2_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 


 




pdf( file = paste(myOutDir_2two, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_500bp), col = km.res4_2two_sub2_500bp$cluster )

plot( t(mat_2two_sub2_500bp), col = km.res4_2two_sub2_500bp$cluster)
points(km.res4_2two_sub2_500bp$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_500bp), col = km.res4_2two_sub2_500bp$cluster)
points(km.res4_2two_sub2_500bp$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_2two_sub2_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "euclidean")
res.dist2_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "maximum")
res.dist3_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "manhattan")
res.dist4_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "canberra")
res.dist5_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "binary")
res.dist6_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "minkowski")
res.dist7_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "pearson")
res.dist8_2two_sub2_500bp <- get_dist( t(mat_2two_sub2_500bp) ,   method = "spearman")
 



pdf( file = paste(myOutDir_2two, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_2two_sub2_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_2two, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_2two_sub2_500bp ) 
print("################### maximum: ")
print(res.dist2_2two_sub2_500bp ) 
print("################### manhattan: ")
print(res.dist3_2two_sub2_500bp ) 
print("################### canberra: ")
print(res.dist4_2two_sub2_500bp ) 
print("################### binary: ")
print(res.dist5_2two_sub2_500bp ) 
print("################### minkowski: ")
print(res.dist6_2two_sub2_500bp ) 
print("################### pearson: ")
print(res.dist7_2two_sub2_500bp ) 
print("################### spearman: ")
print(res.dist8_2two_sub2_500bp ) 
sink() 


# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "ward.D" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "ward.D" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "ward.D" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "ward.D" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "ward.D" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "ward.D" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "ward.D" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "ward.D" )  

pdf( file = paste(myOutDir_2two, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "ward.D2" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "ward.D2" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "ward.D2" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "ward.D2" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "ward.D2" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "ward.D2" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "ward.D2" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "ward.D2" )  

pdf( file = paste(myOutDir_2two, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "single" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "single" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "single" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "single" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "single" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "single" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "single" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "single" )  

pdf( file = paste(myOutDir_2two, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "complete" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "complete" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "complete" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "complete" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "complete" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "complete" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "complete" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "complete" )  

pdf( file = paste(myOutDir_2two, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "average" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "average" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "average" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "average" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "average" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "average" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "average" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "average" )  

pdf( file = paste(myOutDir_2two, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "mcquitty" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "mcquitty" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "mcquitty" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "mcquitty" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "mcquitty" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "mcquitty" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "mcquitty" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "mcquitty" )  

pdf( file = paste(myOutDir_2two, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "median" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "median" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "median" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "median" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "median" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "median" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "median" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "median" )  

pdf( file = paste(myOutDir_2two, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_2two_sub2_500bp <- hclust(res.dist1_2two_sub2_500bp, method = "centroid" )   
res.hc2_2two_sub2_500bp <- hclust(res.dist2_2two_sub2_500bp, method = "centroid" )   
res.hc3_2two_sub2_500bp <- hclust(res.dist3_2two_sub2_500bp, method = "centroid" )   
res.hc4_2two_sub2_500bp <- hclust(res.dist4_2two_sub2_500bp, method = "centroid" )   
res.hc5_2two_sub2_500bp <- hclust(res.dist5_2two_sub2_500bp, method = "centroid" )   
res.hc6_2two_sub2_500bp <- hclust(res.dist6_2two_sub2_500bp, method = "centroid" )   
res.hc7_2two_sub2_500bp <- hclust(res.dist7_2two_sub2_500bp, method = "centroid" )   
res.hc8_2two_sub2_500bp <- hclust(res.dist8_2two_sub2_500bp, method = "centroid" )  

pdf( file = paste(myOutDir_2two, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp )
fviz_dend(res.hc2_2two_sub2_500bp )
fviz_dend(res.hc3_2two_sub2_500bp )
fviz_dend(res.hc4_2two_sub2_500bp )
fviz_dend(res.hc5_2two_sub2_500bp )
fviz_dend(res.hc6_2two_sub2_500bp )
fviz_dend(res.hc7_2two_sub2_500bp )
fviz_dend(res.hc8_2two_sub2_500bp )
dev.off() 


pdf( file = paste(myOutDir_2two, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_500bp, label_cols = myType3_color )
dev.off() 



##################################################################################################################
##################################################################################################################
















######################################################################################################################################################
######################################################################################################################################################
myOutDir_3three = paste(myOutDir, "/3-3kbBin",  sep="");
if( ! file.exists(myOutDir_3three) ) { dir.create(myOutDir_3three, recursive = TRUE) }


#####
tiles_3three_sub3_3kb = tileMethylCounts( myobj, win.size=3000, step.size=3000) ## 3000bp bin
meth_3three_sub3_3kb  = unite( tiles_3three_sub3_3kb, destrand=FALSE   )   ## 100% overlap
mat_3three_sub3_3kb   = percMethylation( meth_3three_sub3_3kb )


sink( file=paste(myOutDir_3three , "0-dimensions-3kbBin.txt", sep="/")  )
print( tiles_3three_sub3_3kb )
print("#########dimensions:")
print( dim(meth_3three_sub3_3kb)  )   
print( dim(mat_3three_sub3_3kb)   )
sink()

meth_3three_sub3_3kb[1:10,]
mat_3three_sub3_3kb[1:10,]



write.table(meth_3three_sub3_3kb , 
            file = paste(myOutDir_3three,   "1A-meth-3kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub3_3kb , 
            file = paste(myOutDir_3three,   "1B-mat-3kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_3three, "2A-MethylationStats-3kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub3_3kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three, "2B-MethylationStats-3kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub3_3kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_3three, "3A-CoverageStats-3kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub3_3kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three, "3B-CoverageStats-3kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub3_3kb[[i]] )  )
}
sink()








pdf( file=paste(myOutDir_3three, "4A-clusterSamples-default-3kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="ward",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="ward",     plot=TRUE)

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="single",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="single",     plot=TRUE)

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="complete",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="complete",     plot=TRUE)

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="average",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="average",     plot=TRUE)

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="mcquitty",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="mcquitty",     plot=TRUE)

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="median",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="median",     plot=TRUE)

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="centroid",     plot=TRUE)
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()








pdf( file=paste(myOutDir_3three, "4B-clusterSamples-noFilter-3kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()






pdf( file=paste(myOutDir_3three, "4C-clusterSamples-sd90perc-3kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()





pdf( file=paste(myOutDir_3three, "4D-clusterSamples-sd99perc-3kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()








pdf( file=paste(myOutDir_3three, "4E-clusterSamples-sd80perc-3kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_3three_sub3_3kb, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_3three_sub3_3kb, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()





pdf( file=paste(myOutDir_3three, "5-PCA-3kbBin.pdf", sep="/") , width=15, height=5  )
PCASamples(meth_3three_sub3_3kb , screeplot=TRUE)
PCASamples(meth_3three_sub3_3kb )
dev.off()

#####################





PCA_3three_sub3_3kb <- prcomp( t(mat_3three_sub3_3kb)  )
names(PCA_3three_sub3_3kb)

sink( file = paste(myOutDir_3three,  "6A-PCA-3kbBin.txt",  sep="/") )
print(PCA_3three_sub3_3kb)
sink()

sink( file = paste(myOutDir_3three,  "6B-PCA-summary-3kbBin.txt",  sep="/") )
summary(PCA_3three_sub3_3kb)
sink()


sink( file = paste(myOutDir_3three,  "6C-PCA-all-3kbBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_3three_sub3_3kb$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_3three_sub3_3kb$rotation)
print("####################### myPCA2$center #########################")
print(PCA_3three_sub3_3kb$center)
print("####################### myPCA2$scale #########################")
print(PCA_3three_sub3_3kb$scale)
print("####################### myPCA2$x #########################")
print(PCA_3three_sub3_3kb$x)
sink()


pdf( file=paste(myOutDir_3three,   "6D-PCA-info-3kbBin.pdf", sep="/")  )
plot(PCA_3three_sub3_3kb, type="lines")
fviz_eig(PCA_3three_sub3_3kb)
dev.off() 




my_fviz_pca_ind1_3three_sub3_3kb <- fviz_pca_ind(PCA_3three_sub3_3kb,
                                                   col.ind = "cos2", # Color by the quality of representation
                                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                   repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_3three_sub3_3kb <- fviz_pca_ind(PCA_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   addEllipses = TRUE, # Concentration ellipses
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE
)
my_fviz_pca_ind3_3three_sub3_3kb <- fviz_pca_ind(PCA_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ) , # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   addEllipses = TRUE, # Concentration ellipses
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE, 
                                                   label = "none", 
                                                   alpha.ind = 1
)
my_fviz_pca_ind4_3three_sub3_3kb <- fviz_pca_ind(PCA_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   #legend.title = "Groups",
                                                   repel = TRUE, 
                                                   label = "none", 
                                                   alpha.ind = 1
)


svg(file=paste(myOutDir_3three, "7A-PCA-2D-1-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "7B-PCA-2D-2-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "7C-PCA-2D-3-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "7D-PCA-2D-4-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_3three_sub3_3kb)
dev.off() 



#############################


PCA_3three_sub3_3kb_matrix <- PCA_3three_sub3_3kb$x
dim( PCA_3three_sub3_3kb_matrix )

PCA_3three_sub3_3kb_Contri  <- (PCA_3three_sub3_3kb$sdev)^2
PCA_3three_sub3_3kb_Contri  <- PCA_3three_sub3_3kb_Contri/sum(PCA_3three_sub3_3kb_Contri)
PCA_3three_sub3_3kb_Contri  <- PCA_3three_sub3_3kb_Contri * 100
PCA_3three_sub3_3kb_Contri  <- round(PCA_3three_sub3_3kb_Contri, 2)

label1_3three_sub3_3kb <-   paste( "PC1 ",  "(", PCA_3three_sub3_3kb_Contri[1], "%)", sep="" )
label2_3three_sub3_3kb <-   paste( "PC2 ",  "(", PCA_3three_sub3_3kb_Contri[2], "%)", sep="" )
label3_3three_sub3_3kb <-   paste( "PC3 ",  "(", PCA_3three_sub3_3kb_Contri[3], "%)", sep="" )
label1_3three_sub3_3kb  
label2_3three_sub3_3kb  
label3_3three_sub3_3kb  


myLabel = myType1
dataframeA_3three_sub3_3kb  <- data.frame( as.data.frame(PCA_3three_sub3_3kb_matrix), myType2, myType3, myLabel   ) 
dataframeA_3three_sub3_3kb 


FigureTemp1_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="8A-PCA-PC1-PC2-3kbBin",  height1=3,  width1=4.8)


FigureTemp2_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="8B-PCA-PC1-PC2-alpha-3kbBin",   height1=3,  width1=4.8)


FigureTemp3_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="8C-PCA-PC1-PC2-smallDot-3kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="8D-PCA-PC1-PC2-big-3kbBin",   height1=3,  width1=4.8)



FigureTemp5_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="8E-PCA-PC1-PC2-text-3kbBin",   height1=3,  width1=4.8)


FigureTemp6_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="8F-PCA-PC1-PC2-text2-3kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="9A-PCA-PC1-PC3-3kbBin",   height1=3,  width1=4.8)


FigureTemp12_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="9B-PCA-PC1-PC3-alpha-3kbBin",   height1=3,  width1=4.8)


FigureTemp13_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="9C-PCA-PC1-PC3-smallDot-3kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="9D-PCA-PC1-PC3-big-3kbBin",  height1=3,  width1=4.8)



FigureTemp15_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="9E-PCA-PC1-PC3-text-3kbBin",   height1=3,  width1=4.8)


FigureTemp16_3three_sub3_3kb  <- ggplot( data = dataframeA_3three_sub3_3kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="9F-PCA-PC1-PC3-text2-3kbBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_3three, "10A_PCA-3d-3kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_3three_sub3_3kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_3kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_3three_sub3_3kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_3kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_3three_sub3_3kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_3kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_3three, "10B_PCA-3d-3kbBin-withLabels.pdf",  sep="/") )

s3d1_3three_sub3_3kb <- scatterplot3d( 
  dataframeA_3three_sub3_3kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_3kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_3three_sub3_3kb$xyz.convert(dataframeA_3three_sub3_3kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_3three_sub3_3kb <- scatterplot3d( 
  dataframeA_3three_sub3_3kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_3kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_3three_sub3_3kb$xyz.convert(dataframeA_3three_sub3_3kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_3three_sub3_3kb <- scatterplot3d( 
  dataframeA_3three_sub3_3kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_3kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_3three_sub3_3kb$xyz.convert(dataframeA_3three_sub3_3kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "11A_PCA-3d-3kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "11B_PCA-3d-3kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_3three, "12_PCA-3d-3kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )



##
scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )



##
scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeA_3three_sub3_3kb[,1], y = dataframeA_3three_sub3_3kb[,2], z = dataframeA_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.


myres.pca_3three_sub3_3kb  <- PCA(t(mat_3three_sub3_3kb),   graph = FALSE)

sink( file = paste(myOutDir_3three, "20_PCA-info-3kbBin-byPCA.txt",  sep="/") )
print( myres.pca_3three_sub3_3kb )
print( "#################################" )
print( summary(myres.pca_3three_sub3_3kb) )
print( "#################################"  )
myeig.val_3three_sub3_3kb <- get_eigenvalue( myres.pca_3three_sub3_3kb )
myeig.val_3three_sub3_3kb
sink() 


pdf( file = paste(myOutDir_3three, "21_PCA-screePlot-3kbBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_3three_sub3_3kb, addlabels = TRUE )
fviz_screeplot(X=myres.pca_3three_sub3_3kb, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_3three_sub3_3kb, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_3three_sub3_3kb, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_3three_sub3_3kb, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_3three_sub3_3kb <- fviz_pca_ind(myres.pca_3three_sub3_3kb,
                                                   col.ind = "cos2", # Color by the quality of representation
                                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                   repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_3three_sub3_3kb <- fviz_pca_ind(myres.pca_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   addEllipses = TRUE, # Concentration ellipses
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE
)
my_fviz_pca_ind3_3three_sub3_3kb <- fviz_pca_ind(myres.pca_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ) , # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   addEllipses = TRUE, # Concentration ellipses
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE, 
                                                   label = "none", 
                                                   alpha.ind = 1
)
my_fviz_pca_ind4_3three_sub3_3kb <- fviz_pca_ind(myres.pca_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   #legend.title = "Groups",
                                                   repel = TRUE, 
                                                   label = "none", 
                                                   alpha.ind = 1
)
my_fviz_pca_ind5_3three_sub3_3kb <- fviz_pca_ind(myres.pca_3three_sub3_3kb,
                                                   col.ind =  as.factor( myTreatment ), # color by groups
                                                   #palette = c("#00AFBB",  "#FC4E07"),
                                                   ellipse.type = "confidence",
                                                   legend.title = "Groups",
                                                   repel = TRUE
)

svg(file=paste(myOutDir_3three, "22A-PCA-2D-1-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22B-PCA-2D-2-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22C-PCA-2D-3-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22D-PCA-2D-4-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_3three_sub3_3kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22E-PCA-2D-4-3kbBin.svg", sep="/") )
print(my_fviz_pca_ind5_3three_sub3_3kb)
dev.off() 



#############################


myres.pca_3three_sub3_3kb_matrix <- myres.pca_3three_sub3_3kb$ind$coord
dim( myres.pca_3three_sub3_3kb_matrix )

myres.pca_3three_sub3_3kb_Contri  <- (myres.pca_3three_sub3_3kb$eig)[,2] 
myres.pca_3three_sub3_3kb_Contri  <- myres.pca_3three_sub3_3kb_Contri/sum(myres.pca_3three_sub3_3kb_Contri)
myres.pca_3three_sub3_3kb_Contri  <- myres.pca_3three_sub3_3kb_Contri * 100
myres.pca_3three_sub3_3kb_Contri  <- round(myres.pca_3three_sub3_3kb_Contri, 2)

label1_3three_sub3_3kb <-   paste( "PC1 ",  "(", myres.pca_3three_sub3_3kb_Contri[1], "%)", sep="" )
label2_3three_sub3_3kb <-   paste( "PC2 ",  "(", myres.pca_3three_sub3_3kb_Contri[2], "%)", sep="" )
label3_3three_sub3_3kb <-   paste( "PC3 ",  "(", myres.pca_3three_sub3_3kb_Contri[3], "%)", sep="" )
label1_3three_sub3_3kb  
label2_3three_sub3_3kb  
label3_3three_sub3_3kb  


myLabel = myType1
dataframeB_3three_sub3_3kb  <- data.frame( as.data.frame(myres.pca_3three_sub3_3kb_matrix), myType2, myType3, myLabel   ) 
dataframeB_3three_sub3_3kb 


FigureTemp1_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="23A-PCA-PC1-PC2-3kbBin",  height1=3,  width1=4.8)


FigureTemp2_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="23B-PCA-PC1-PC2-alpha-3kbBin",   height1=3,  width1=4.8)


FigureTemp3_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="23C-PCA-PC1-PC2-smallDot-3kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="23D-PCA-PC1-PC2-big-3kbBin",   height1=3,  width1=4.8)



FigureTemp5_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="23E-PCA-PC1-PC2-text-3kbBin",   height1=3,  width1=4.8)


FigureTemp6_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label2_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="23F-PCA-PC1-PC2-text2-3kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="24A-PCA-PC1-PC3-3kbBin",   height1=3,  width1=4.8)


FigureTemp12_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="24B-PCA-PC1-PC3-alpha-3kbBin",   height1=3,  width1=4.8)


FigureTemp13_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="24C-PCA-PC1-PC3-smallDot-3kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="24D-PCA-PC1-PC3-big-3kbBin",  height1=3,  width1=4.8)



FigureTemp15_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="24E-PCA-PC1-PC3-text-3kbBin",   height1=3,  width1=4.8)


FigureTemp16_3three_sub3_3kb  <- ggplot( data = dataframeB_3three_sub3_3kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_3kb) +   ylab(label3_3three_sub3_3kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_3three_sub3_3kb ,  path1=myOutDir_3three, fileName1="24F-PCA-PC1-PC3-text2-3kbBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_3three, "25A_PCA-3d-3kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_3three_sub3_3kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_3kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_3three_sub3_3kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_3kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_3three_sub3_3kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_3kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_3three, "25B_PCA-3d-3kbBin-withLabels.pdf",  sep="/") )

s3d1_3three_sub3_3kb <- scatterplot3d( 
  dataframeB_3three_sub3_3kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_3kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_3three_sub3_3kb$xyz.convert(dataframeB_3three_sub3_3kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_3three_sub3_3kb <- scatterplot3d( 
  dataframeB_3three_sub3_3kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_3kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_3three_sub3_3kb$xyz.convert(dataframeB_3three_sub3_3kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_3three_sub3_3kb <- scatterplot3d( 
  dataframeB_3three_sub3_3kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_3kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_3three_sub3_3kb$xyz.convert(dataframeB_3three_sub3_3kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "26A_PCA-3d-3kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "26B_PCA-3d-3kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_3three, "27_PCA-3d-3kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )



##
scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )



##
scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )


scatter3D(x = dataframeB_3three_sub3_3kb[,1], y = dataframeB_3three_sub3_3kb[,2], z = dataframeB_3three_sub3_3kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_3kb,  ylab = label2_3three_sub3_3kb,  zlab = label3_3three_sub3_3kb  )

dev.off()






#################



km.res1_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=2, iter.max = 20 )
km.res2_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=3, iter.max = 20  )
km.res3_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=4, iter.max = 20  )
km.res4_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=5, iter.max = 20 )
km.res5_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=6, iter.max = 20 )
km.res6_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=7, iter.max = 20 )
km.res7_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=8, iter.max = 20 )
km.res8_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=9, iter.max = 20 )
km.res9_3three_sub3_3kb  <- kmeans( t(mat_3three_sub3_3kb), centers=10, iter.max = 20 )
km.res10_3three_sub3_3kb <- kmeans( t(mat_3three_sub3_3kb), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_3three, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_3three_sub3_3kb$cluster
print("#################### 3 classes:")
km.res2_3three_sub3_3kb$cluster
print("#################### 4 classes:")
km.res3_3three_sub3_3kb$cluster
print("#################### 5 classes:")
km.res4_3three_sub3_3kb$cluster
print("#################### 6 classes:")
km.res5_3three_sub3_3kb$cluster
print("#################### 7 classes:")
km.res6_3three_sub3_3kb$cluster
print("#################### 8 classes:")
km.res7_3three_sub3_3kb$cluster
print("#################### 9 classes:")
km.res8_3three_sub3_3kb$cluster
print("#################### 10 classes:")
km.res9_3three_sub3_3kb$cluster
print("#################### 11 classes:")
km.res10_3three_sub3_3kb$cluster
sink()



pdf( file = paste(myOutDir_3three, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_3kb), col = km.res1_3three_sub3_3kb$cluster )

plot( t(mat_3three_sub3_3kb), col = km.res1_3three_sub3_3kb$cluster)
points(km.res1_3three_sub3_3kb$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_3kb), col = km.res1_3three_sub3_3kb$cluster)
points(km.res1_3three_sub3_3kb$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_3three_sub3_3kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 






pdf( file = paste(myOutDir_3three, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_3kb), col = km.res2_3three_sub3_3kb$cluster )

plot( t(mat_3three_sub3_3kb), col = km.res2_3three_sub3_3kb$cluster)
points(km.res2_3three_sub3_3kb$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_3kb), col = km.res2_3three_sub3_3kb$cluster)
points(km.res2_3three_sub3_3kb$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_3three_sub3_3kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_3three, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_3kb), col = km.res3_3three_sub3_3kb$cluster )

plot( t(mat_3three_sub3_3kb), col = km.res3_3three_sub3_3kb$cluster)
points(km.res3_3three_sub3_3kb$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_3kb), col = km.res3_3three_sub3_3kb$cluster)
points(km.res3_3three_sub3_3kb$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_3three_sub3_3kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 







pdf( file = paste(myOutDir_3three, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_3kb), col = km.res4_3three_sub3_3kb$cluster )

plot( t(mat_3three_sub3_3kb), col = km.res4_3three_sub3_3kb$cluster)
points(km.res4_3three_sub3_3kb$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_3kb), col = km.res4_3three_sub3_3kb$cluster)
points(km.res4_3three_sub3_3kb$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_3three_sub3_3kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "euclidean")
res.dist2_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "maximum")
res.dist3_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "manhattan")
res.dist4_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "canberra")
res.dist5_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "binary")
res.dist6_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "minkowski")
res.dist7_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "pearson")
res.dist8_3three_sub3_3kb <- get_dist( t(mat_3three_sub3_3kb) ,   method = "spearman")




pdf( file = paste(myOutDir_3three, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_3three_sub3_3kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_3three, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_3three_sub3_3kb ) 
print("################### maximum: ")
print(res.dist2_3three_sub3_3kb ) 
print("################### manhattan: ")
print(res.dist3_3three_sub3_3kb ) 
print("################### canberra: ")
print(res.dist4_3three_sub3_3kb ) 
print("################### binary: ")
print(res.dist5_3three_sub3_3kb ) 
print("################### minkowski: ")
print(res.dist6_3three_sub3_3kb ) 
print("################### pearson: ")
print(res.dist7_3three_sub3_3kb ) 
print("################### spearman: ")
print(res.dist8_3three_sub3_3kb ) 
sink() 


# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "ward.D" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "ward.D" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "ward.D" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "ward.D" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "ward.D" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "ward.D" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "ward.D" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "ward.D" )  

pdf( file = paste(myOutDir_3three, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "ward.D2" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "ward.D2" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "ward.D2" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "ward.D2" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "ward.D2" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "ward.D2" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "ward.D2" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "ward.D2" )  

pdf( file = paste(myOutDir_3three, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "single" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "single" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "single" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "single" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "single" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "single" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "single" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "single" )  

pdf( file = paste(myOutDir_3three, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "complete" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "complete" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "complete" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "complete" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "complete" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "complete" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "complete" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "complete" )  

pdf( file = paste(myOutDir_3three, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "average" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "average" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "average" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "average" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "average" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "average" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "average" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "average" )  

pdf( file = paste(myOutDir_3three, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "mcquitty" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "mcquitty" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "mcquitty" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "mcquitty" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "mcquitty" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "mcquitty" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "mcquitty" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "mcquitty" )  

pdf( file = paste(myOutDir_3three, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "median" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "median" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "median" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "median" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "median" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "median" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "median" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "median" )  

pdf( file = paste(myOutDir_3three, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_3three_sub3_3kb <- hclust(res.dist1_3three_sub3_3kb, method = "centroid" )   
res.hc2_3three_sub3_3kb <- hclust(res.dist2_3three_sub3_3kb, method = "centroid" )   
res.hc3_3three_sub3_3kb <- hclust(res.dist3_3three_sub3_3kb, method = "centroid" )   
res.hc4_3three_sub3_3kb <- hclust(res.dist4_3three_sub3_3kb, method = "centroid" )   
res.hc5_3three_sub3_3kb <- hclust(res.dist5_3three_sub3_3kb, method = "centroid" )   
res.hc6_3three_sub3_3kb <- hclust(res.dist6_3three_sub3_3kb, method = "centroid" )   
res.hc7_3three_sub3_3kb <- hclust(res.dist7_3three_sub3_3kb, method = "centroid" )   
res.hc8_3three_sub3_3kb <- hclust(res.dist8_3three_sub3_3kb, method = "centroid" )  

pdf( file = paste(myOutDir_3three, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_3kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb )
fviz_dend(res.hc2_3three_sub3_3kb )
fviz_dend(res.hc3_3three_sub3_3kb )
fviz_dend(res.hc4_3three_sub3_3kb )
fviz_dend(res.hc5_3three_sub3_3kb )
fviz_dend(res.hc6_3three_sub3_3kb )
fviz_dend(res.hc7_3three_sub3_3kb )
fviz_dend(res.hc8_3three_sub3_3kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_3kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_3kb, label_cols = myType3_color )
dev.off() 



##################################################################################################################
##################################################################################################################






























######################################################################################################################################################
######################################################################################################################################################
myOutDir_4four = paste(myOutDir, "/4-5kbBin",  sep="");
if( ! file.exists(myOutDir_4four) ) { dir.create(myOutDir_4four, recursive = TRUE) }


#####
tiles_4four_sub4_5kb = tileMethylCounts( myobj, win.size=5000, step.size=5000) ## 5000bp bin
meth_4four_sub4_5kb  = unite( tiles_4four_sub4_5kb, destrand=FALSE   )   ## 100% overlap
mat_4four_sub4_5kb   = percMethylation( meth_4four_sub4_5kb )


sink( file=paste(myOutDir_4four , "0-dimensions-5kbBin.txt", sep="/")  )
print( tiles_4four_sub4_5kb )
print("#########dimensions:")
print( dim(meth_4four_sub4_5kb)  )   
print( dim(mat_4four_sub4_5kb)   )
sink()

meth_4four_sub4_5kb[1:10,]
mat_4four_sub4_5kb[1:10,]



write.table(meth_4four_sub4_5kb , 
            file = paste(myOutDir_4four,   "1A-meth-5kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub4_5kb , 
            file = paste(myOutDir_4four,   "1B-mat-5kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")



pdf( file=paste(myOutDir_4four, "2A-MethylationStats-5kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub4_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four, "2B-MethylationStats-5kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub4_5kb[[i]] )  )
}
sink()


pdf( file=paste(myOutDir_4four, "3A-CoverageStats-5kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub4_5kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four, "3B-CoverageStats-5kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub4_5kb[[i]] )  )
}
sink()








pdf( file=paste(myOutDir_4four, "4A-clusterSamples-default-5kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="ward",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="ward",     plot=TRUE)

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="single",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="single",     plot=TRUE)

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="complete",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="complete",     plot=TRUE)

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="average",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="average",     plot=TRUE)

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="mcquitty",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="mcquitty",     plot=TRUE)

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="median",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="median",     plot=TRUE)

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="centroid",     plot=TRUE)
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="centroid",     plot=TRUE)

dev.off()








pdf( file=paste(myOutDir_4four, "4B-clusterSamples-noFilter-5kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="ward",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="single",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="complete",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="average",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="mcquitty",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="median",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="centroid",     sd.filter=FALSE,   sd.threshold=0,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()






pdf( file=paste(myOutDir_4four, "4C-clusterSamples-sd90perc-5kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.9,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()





pdf( file=paste(myOutDir_4four, "4D-clusterSamples-sd99perc-5kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.99,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()








pdf( file=paste(myOutDir_4four, "4E-clusterSamples-sd80perc-5kbBin.pdf", sep="/") , width=15, height=5  )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="ward",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="single",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="complete",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="average",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="mcquitty",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="median",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

clusterSamples(meth_4four_sub4_5kb, dist="correlation", method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="euclidean",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="maximum",     method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="manhattan",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="canberra",    method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="binary",      method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )
clusterSamples(meth_4four_sub4_5kb, dist="minkowski",   method="centroid",     sd.filter=TRUE,   sd.threshold=0.8,                filterByQuantile=TRUE,   plot=TRUE )

dev.off()





pdf( file=paste(myOutDir_4four, "5-PCA-5kbBin.pdf", sep="/") , width=15, height=5  )
PCASamples(meth_4four_sub4_5kb , screeplot=TRUE)
PCASamples(meth_4four_sub4_5kb )
dev.off()

#####################





PCA_4four_sub4_5kb <- prcomp( t(mat_4four_sub4_5kb)  )
names(PCA_4four_sub4_5kb)

sink( file = paste(myOutDir_4four,  "6A-PCA-5kbBin.txt",  sep="/") )
print(PCA_4four_sub4_5kb)
sink()

sink( file = paste(myOutDir_4four,  "6B-PCA-summary-5kbBin.txt",  sep="/") )
summary(PCA_4four_sub4_5kb)
sink()


sink( file = paste(myOutDir_4four,  "6C-PCA-all-5kbBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_4four_sub4_5kb$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_4four_sub4_5kb$rotation)
print("####################### myPCA2$center #########################")
print(PCA_4four_sub4_5kb$center)
print("####################### myPCA2$scale #########################")
print(PCA_4four_sub4_5kb$scale)
print("####################### myPCA2$x #########################")
print(PCA_4four_sub4_5kb$x)
sink()


pdf( file=paste(myOutDir_4four,   "6D-PCA-info-5kbBin.pdf", sep="/")  )
plot(PCA_4four_sub4_5kb, type="lines")
fviz_eig(PCA_4four_sub4_5kb)
dev.off() 




my_fviz_pca_ind1_4four_sub4_5kb <- fviz_pca_ind(PCA_4four_sub4_5kb,
                                                col.ind = "cos2", # Color by the quality of representation
                                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_4four_sub4_5kb <- fviz_pca_ind(PCA_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                addEllipses = TRUE, # Concentration ellipses
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE
)
my_fviz_pca_ind3_4four_sub4_5kb <- fviz_pca_ind(PCA_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ) , # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                addEllipses = TRUE, # Concentration ellipses
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE, 
                                                label = "none", 
                                                alpha.ind = 1
)
my_fviz_pca_ind4_4four_sub4_5kb <- fviz_pca_ind(PCA_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                #legend.title = "Groups",
                                                repel = TRUE, 
                                                label = "none", 
                                                alpha.ind = 1
)


svg(file=paste(myOutDir_4four, "7A-PCA-2D-1-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "7B-PCA-2D-2-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "7C-PCA-2D-3-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "7D-PCA-2D-4-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_4four_sub4_5kb)
dev.off() 



#############################


PCA_4four_sub4_5kb_matrix <- PCA_4four_sub4_5kb$x
dim( PCA_4four_sub4_5kb_matrix )

PCA_4four_sub4_5kb_Contri  <- (PCA_4four_sub4_5kb$sdev)^2
PCA_4four_sub4_5kb_Contri  <- PCA_4four_sub4_5kb_Contri/sum(PCA_4four_sub4_5kb_Contri)
PCA_4four_sub4_5kb_Contri  <- PCA_4four_sub4_5kb_Contri * 100
PCA_4four_sub4_5kb_Contri  <- round(PCA_4four_sub4_5kb_Contri, 2)

label1_4four_sub4_5kb <-   paste( "PC1 ",  "(", PCA_4four_sub4_5kb_Contri[1], "%)", sep="" )
label2_4four_sub4_5kb <-   paste( "PC2 ",  "(", PCA_4four_sub4_5kb_Contri[2], "%)", sep="" )
label3_4four_sub4_5kb <-   paste( "PC3 ",  "(", PCA_4four_sub4_5kb_Contri[3], "%)", sep="" )
label1_4four_sub4_5kb  
label2_4four_sub4_5kb  
label3_4four_sub4_5kb  


myLabel = myType1
dataframeA_4four_sub4_5kb  <- data.frame( as.data.frame(PCA_4four_sub4_5kb_matrix), myType2, myType3, myLabel   ) 
dataframeA_4four_sub4_5kb 


FigureTemp1_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="8A-PCA-PC1-PC2-5kbBin",  height1=3,  width1=4.8)


FigureTemp2_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="8B-PCA-PC1-PC2-alpha-5kbBin",   height1=3,  width1=4.8)


FigureTemp3_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="8C-PCA-PC1-PC2-smallDot-5kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="8D-PCA-PC1-PC2-big-5kbBin",   height1=3,  width1=4.8)



FigureTemp5_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="8E-PCA-PC1-PC2-text-5kbBin",   height1=3,  width1=4.8)


FigureTemp6_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="8F-PCA-PC1-PC2-text2-5kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="9A-PCA-PC1-PC3-5kbBin",   height1=3,  width1=4.8)


FigureTemp12_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="9B-PCA-PC1-PC3-alpha-5kbBin",   height1=3,  width1=4.8)


FigureTemp13_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="9C-PCA-PC1-PC3-smallDot-5kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="9D-PCA-PC1-PC3-big-5kbBin",  height1=3,  width1=4.8)



FigureTemp15_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="9E-PCA-PC1-PC3-text-5kbBin",   height1=3,  width1=4.8)


FigureTemp16_4four_sub4_5kb  <- ggplot( data = dataframeA_4four_sub4_5kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="9F-PCA-PC1-PC3-text2-5kbBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_4four, "10A_PCA-3d-5kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_4four_sub4_5kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_5kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_4four_sub4_5kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_5kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_4four_sub4_5kb[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_5kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_4four, "10B_PCA-3d-5kbBin-withLabels.pdf",  sep="/") )

s3d1_4four_sub4_5kb <- scatterplot3d( 
  dataframeA_4four_sub4_5kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_5kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_4four_sub4_5kb$xyz.convert(dataframeA_4four_sub4_5kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_4four_sub4_5kb <- scatterplot3d( 
  dataframeA_4four_sub4_5kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_5kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_4four_sub4_5kb$xyz.convert(dataframeA_4four_sub4_5kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_4four_sub4_5kb <- scatterplot3d( 
  dataframeA_4four_sub4_5kb[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_5kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_4four_sub4_5kb$xyz.convert(dataframeA_4four_sub4_5kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "11A_PCA-3d-5kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "11B_PCA-3d-5kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_4four, "12_PCA-3d-5kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )



##
scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )



##
scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeA_4four_sub4_5kb[,1], y = dataframeA_4four_sub4_5kb[,2], z = dataframeA_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.


myres.pca_4four_sub4_5kb  <- PCA(t(mat_4four_sub4_5kb),   graph = FALSE)

sink( file = paste(myOutDir_4four, "20_PCA-info-5kbBin-byPCA.txt",  sep="/") )
print( myres.pca_4four_sub4_5kb )
print( "#################################" )
print( summary(myres.pca_4four_sub4_5kb) )
print( "#################################"  )
myeig.val_4four_sub4_5kb <- get_eigenvalue( myres.pca_4four_sub4_5kb )
myeig.val_4four_sub4_5kb
sink() 


pdf( file = paste(myOutDir_4four, "21_PCA-screePlot-5kbBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_4four_sub4_5kb, addlabels = TRUE )
fviz_screeplot(X=myres.pca_4four_sub4_5kb, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_4four_sub4_5kb, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_4four_sub4_5kb, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_4four_sub4_5kb, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_4four_sub4_5kb <- fviz_pca_ind(myres.pca_4four_sub4_5kb,
                                                col.ind = "cos2", # Color by the quality of representation
                                                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_4four_sub4_5kb <- fviz_pca_ind(myres.pca_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                addEllipses = TRUE, # Concentration ellipses
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE
)
my_fviz_pca_ind3_4four_sub4_5kb <- fviz_pca_ind(myres.pca_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ) , # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                addEllipses = TRUE, # Concentration ellipses
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE, 
                                                label = "none", 
                                                alpha.ind = 1
)
my_fviz_pca_ind4_4four_sub4_5kb <- fviz_pca_ind(myres.pca_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                #legend.title = "Groups",
                                                repel = TRUE, 
                                                label = "none", 
                                                alpha.ind = 1
)
my_fviz_pca_ind5_4four_sub4_5kb <- fviz_pca_ind(myres.pca_4four_sub4_5kb,
                                                col.ind =  as.factor( myTreatment ), # color by groups
                                                #palette = c("#00AFBB",  "#FC4E07"),
                                                ellipse.type = "confidence",
                                                legend.title = "Groups",
                                                repel = TRUE
)

svg(file=paste(myOutDir_4four, "22A-PCA-2D-1-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "22B-PCA-2D-2-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "22C-PCA-2D-3-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "22D-PCA-2D-4-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_4four_sub4_5kb)
dev.off() 

svg(file=paste(myOutDir_4four, "22E-PCA-2D-4-5kbBin.svg", sep="/") )
print(my_fviz_pca_ind5_4four_sub4_5kb)
dev.off() 



#############################


myres.pca_4four_sub4_5kb_matrix <- myres.pca_4four_sub4_5kb$ind$coord
dim( myres.pca_4four_sub4_5kb_matrix )

myres.pca_4four_sub4_5kb_Contri  <- (myres.pca_4four_sub4_5kb$eig)[,2] 
myres.pca_4four_sub4_5kb_Contri  <- myres.pca_4four_sub4_5kb_Contri/sum(myres.pca_4four_sub4_5kb_Contri)
myres.pca_4four_sub4_5kb_Contri  <- myres.pca_4four_sub4_5kb_Contri * 100
myres.pca_4four_sub4_5kb_Contri  <- round(myres.pca_4four_sub4_5kb_Contri, 2)

label1_4four_sub4_5kb <-   paste( "PC1 ",  "(", myres.pca_4four_sub4_5kb_Contri[1], "%)", sep="" )
label2_4four_sub4_5kb <-   paste( "PC2 ",  "(", myres.pca_4four_sub4_5kb_Contri[2], "%)", sep="" )
label3_4four_sub4_5kb <-   paste( "PC3 ",  "(", myres.pca_4four_sub4_5kb_Contri[3], "%)", sep="" )
label1_4four_sub4_5kb  
label2_4four_sub4_5kb  
label3_4four_sub4_5kb  


myLabel = myType1
dataframeB_4four_sub4_5kb  <- data.frame( as.data.frame(myres.pca_4four_sub4_5kb_matrix), myType2, myType3, myLabel   ) 
dataframeB_4four_sub4_5kb 


FigureTemp1_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="23A-PCA-PC1-PC2-5kbBin",  height1=3,  width1=4.8)


FigureTemp2_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="23B-PCA-PC1-PC2-alpha-5kbBin",   height1=3,  width1=4.8)


FigureTemp3_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="23C-PCA-PC1-PC2-smallDot-5kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="23D-PCA-PC1-PC2-big-5kbBin",   height1=3,  width1=4.8)



FigureTemp5_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="23E-PCA-PC1-PC2-text-5kbBin",   height1=3,  width1=4.8)


FigureTemp6_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label2_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="23F-PCA-PC1-PC2-text2-5kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="24A-PCA-PC1-PC3-5kbBin",   height1=3,  width1=4.8)


FigureTemp12_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="24B-PCA-PC1-PC3-alpha-5kbBin",   height1=3,  width1=4.8)


FigureTemp13_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="24C-PCA-PC1-PC3-smallDot-5kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="24D-PCA-PC1-PC3-big-5kbBin",  height1=3,  width1=4.8)



FigureTemp15_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="24E-PCA-PC1-PC3-text-5kbBin",   height1=3,  width1=4.8)


FigureTemp16_4four_sub4_5kb  <- ggplot( data = dataframeB_4four_sub4_5kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_5kb) +   ylab(label3_4four_sub4_5kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_4four_sub4_5kb ,  path1=myOutDir_4four, fileName1="24F-PCA-PC1-PC3-text2-5kbBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_4four, "25A_PCA-3d-5kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_4four_sub4_5kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_5kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_4four_sub4_5kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_5kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_4four_sub4_5kb[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_5kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_4four, "25B_PCA-3d-5kbBin-withLabels.pdf",  sep="/") )

s3d1_4four_sub4_5kb <- scatterplot3d( 
  dataframeB_4four_sub4_5kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_5kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_4four_sub4_5kb$xyz.convert(dataframeB_4four_sub4_5kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_4four_sub4_5kb <- scatterplot3d( 
  dataframeB_4four_sub4_5kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_5kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_4four_sub4_5kb$xyz.convert(dataframeB_4four_sub4_5kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_4four_sub4_5kb <- scatterplot3d( 
  dataframeB_4four_sub4_5kb[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_5kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_4four_sub4_5kb$xyz.convert(dataframeB_4four_sub4_5kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "26A_PCA-3d-5kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "26B_PCA-3d-5kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_4four, "27_PCA-3d-5kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )



##
scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )



##
scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )


scatter3D(x = dataframeB_4four_sub4_5kb[,1], y = dataframeB_4four_sub4_5kb[,2], z = dataframeB_4four_sub4_5kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_5kb,  ylab = label2_4four_sub4_5kb,  zlab = label3_4four_sub4_5kb  )

dev.off()






#################



km.res1_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=2, iter.max = 20 )
km.res2_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=3, iter.max = 20  )
km.res3_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=4, iter.max = 20  )
km.res4_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=5, iter.max = 20 )
km.res5_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=6, iter.max = 20 )
km.res6_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=7, iter.max = 20 )
km.res7_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=8, iter.max = 20 )
km.res8_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=9, iter.max = 20 )
km.res9_4four_sub4_5kb  <- kmeans( t(mat_4four_sub4_5kb), centers=10, iter.max = 20 )
km.res10_4four_sub4_5kb <- kmeans( t(mat_4four_sub4_5kb), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_4four, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_4four_sub4_5kb$cluster
print("#################### 3 classes:")
km.res2_4four_sub4_5kb$cluster
print("#################### 4 classes:")
km.res3_4four_sub4_5kb$cluster
print("#################### 5 classes:")
km.res4_4four_sub4_5kb$cluster
print("#################### 6 classes:")
km.res5_4four_sub4_5kb$cluster
print("#################### 7 classes:")
km.res6_4four_sub4_5kb$cluster
print("#################### 8 classes:")
km.res7_4four_sub4_5kb$cluster
print("#################### 9 classes:")
km.res8_4four_sub4_5kb$cluster
print("#################### 10 classes:")
km.res9_4four_sub4_5kb$cluster
print("#################### 11 classes:")
km.res10_4four_sub4_5kb$cluster
sink()



pdf( file = paste(myOutDir_4four, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_5kb), col = km.res1_4four_sub4_5kb$cluster )

plot( t(mat_4four_sub4_5kb), col = km.res1_4four_sub4_5kb$cluster)
points(km.res1_4four_sub4_5kb$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_5kb), col = km.res1_4four_sub4_5kb$cluster)
points(km.res1_4four_sub4_5kb$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_4four_sub4_5kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 






pdf( file = paste(myOutDir_4four, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_5kb), col = km.res2_4four_sub4_5kb$cluster )

plot( t(mat_4four_sub4_5kb), col = km.res2_4four_sub4_5kb$cluster)
points(km.res2_4four_sub4_5kb$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_5kb), col = km.res2_4four_sub4_5kb$cluster)
points(km.res2_4four_sub4_5kb$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_4four_sub4_5kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_4four, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_5kb), col = km.res3_4four_sub4_5kb$cluster )

plot( t(mat_4four_sub4_5kb), col = km.res3_4four_sub4_5kb$cluster)
points(km.res3_4four_sub4_5kb$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_5kb), col = km.res3_4four_sub4_5kb$cluster)
points(km.res3_4four_sub4_5kb$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_4four_sub4_5kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 







pdf( file = paste(myOutDir_4four, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_5kb), col = km.res4_4four_sub4_5kb$cluster )

plot( t(mat_4four_sub4_5kb), col = km.res4_4four_sub4_5kb$cluster)
points(km.res4_4four_sub4_5kb$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_5kb), col = km.res4_4four_sub4_5kb$cluster)
points(km.res4_4four_sub4_5kb$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_4four_sub4_5kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "euclidean")
res.dist2_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "maximum")
res.dist3_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "manhattan")
res.dist4_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "canberra")
res.dist5_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "binary")
res.dist6_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "minkowski")
res.dist7_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "pearson")
res.dist8_4four_sub4_5kb <- get_dist( t(mat_4four_sub4_5kb) ,   method = "spearman")




pdf( file = paste(myOutDir_4four, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_4four_sub4_5kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_4four, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_4four_sub4_5kb ) 
print("################### maximum: ")
print(res.dist2_4four_sub4_5kb ) 
print("################### manhattan: ")
print(res.dist3_4four_sub4_5kb ) 
print("################### canberra: ")
print(res.dist4_4four_sub4_5kb ) 
print("################### binary: ")
print(res.dist5_4four_sub4_5kb ) 
print("################### minkowski: ")
print(res.dist6_4four_sub4_5kb ) 
print("################### pearson: ")
print(res.dist7_4four_sub4_5kb ) 
print("################### spearman: ")
print(res.dist8_4four_sub4_5kb ) 
sink() 


# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "ward.D" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "ward.D" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "ward.D" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "ward.D" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "ward.D" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "ward.D" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "ward.D" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "ward.D" )  

pdf( file = paste(myOutDir_4four, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "ward.D2" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "ward.D2" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "ward.D2" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "ward.D2" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "ward.D2" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "ward.D2" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "ward.D2" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "ward.D2" )  

pdf( file = paste(myOutDir_4four, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "single" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "single" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "single" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "single" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "single" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "single" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "single" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "single" )  

pdf( file = paste(myOutDir_4four, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "complete" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "complete" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "complete" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "complete" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "complete" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "complete" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "complete" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "complete" )  

pdf( file = paste(myOutDir_4four, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "average" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "average" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "average" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "average" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "average" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "average" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "average" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "average" )  

pdf( file = paste(myOutDir_4four, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "mcquitty" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "mcquitty" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "mcquitty" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "mcquitty" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "mcquitty" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "mcquitty" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "mcquitty" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "mcquitty" )  

pdf( file = paste(myOutDir_4four, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "median" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "median" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "median" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "median" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "median" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "median" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "median" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "median" )  

pdf( file = paste(myOutDir_4four, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_4four_sub4_5kb <- hclust(res.dist1_4four_sub4_5kb, method = "centroid" )   
res.hc2_4four_sub4_5kb <- hclust(res.dist2_4four_sub4_5kb, method = "centroid" )   
res.hc3_4four_sub4_5kb <- hclust(res.dist3_4four_sub4_5kb, method = "centroid" )   
res.hc4_4four_sub4_5kb <- hclust(res.dist4_4four_sub4_5kb, method = "centroid" )   
res.hc5_4four_sub4_5kb <- hclust(res.dist5_4four_sub4_5kb, method = "centroid" )   
res.hc6_4four_sub4_5kb <- hclust(res.dist6_4four_sub4_5kb, method = "centroid" )   
res.hc7_4four_sub4_5kb <- hclust(res.dist7_4four_sub4_5kb, method = "centroid" )   
res.hc8_4four_sub4_5kb <- hclust(res.dist8_4four_sub4_5kb, method = "centroid" )  

pdf( file = paste(myOutDir_4four, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_5kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb )
fviz_dend(res.hc2_4four_sub4_5kb )
fviz_dend(res.hc3_4four_sub4_5kb )
fviz_dend(res.hc4_4four_sub4_5kb )
fviz_dend(res.hc5_4four_sub4_5kb )
fviz_dend(res.hc6_4four_sub4_5kb )
fviz_dend(res.hc7_4four_sub4_5kb )
fviz_dend(res.hc8_4four_sub4_5kb )
dev.off() 


pdf( file = paste(myOutDir_4four, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_5kb, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_5kb, label_cols = myType3_color )
dev.off() 



##################################################################################################################
##################################################################################################################






