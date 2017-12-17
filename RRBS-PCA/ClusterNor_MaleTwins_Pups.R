## All samples must be overlapped. (100% overlap)
##
## Example:  
## Rscript  Cluster_MaleTwins_Pups.R     11B_splitXY/100A-rmXY/5_cov50reads     1011B_splitXY/100A-rmXY/5_cov50reads/Cluster_MaleTwins_Pups          
         
args <- commandArgs(TRUE)
print("args: ")
print(args[1])   
print(args[2])     
print("#############")

inputDir   = args[1];     ## the path of input file
myOutDir   = args[2];     ## the path of output file

# inputDir =  "11B_splitXY/100A-rmXY/1_cov10reads"
# myOutDir =  "1011B_splitXY/100A-rmXY/1_cov10reads/Cluster_MaleTwins_Pups"





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





Files_NC <- c(
  paste( inputDir,   "67_E24C-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir,   "67_E24D-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir,   "68_E56C-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir,   "68_E56D-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir,   "69_NC-E123-C-Boy_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir,   "69_NC-E123-D-Boy_Rep1.bismark.cov",  sep="/" ) 
)
mySampleID_NC <- rep( x="NC", times=length(Files_NC) )
for(i in c(1:length(mySampleID_NC)) ) {
  mySampleID_NC[i] = paste(mySampleID_NC[i], i, sep="_")
}
myTreatment_NC <- rep( x=0,  times=length(Files_NC) )
Sex_NC = rep( x="boy",  times=length(Files_NC) )  
for(i in c(1:length(Sex_NC)) ) {
  Sex_NC[i] = "boy"
  if(i>=4) { Sex_NC[i] = "boy" }
}

  

Files_IVF_fresh <- c(
  paste( inputDir,   "70_E113C-boy-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir,   "70_E113D-boy-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir,   "71_ART-W58-C-Boy_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "71_ART-W58-D-Boy_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "72_ART-W779-C-Boy_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir,   "72_W779D-ART-boy_Rep1.bismark.cov",    sep="/" ) 
)
mySampleID_IVF_fresh <- rep( x="IVF_fresh", times=length(Files_IVF_fresh) )
for(i in c(1:length(mySampleID_IVF_fresh)) ) {
  mySampleID_IVF_fresh[i] = paste(mySampleID_IVF_fresh[i], i, sep="_")
}
myTreatment_IVF_fresh <- rep( x=1,  times=length(Files_IVF_fresh) )
Sex_IVF_fresh = rep( x="boy",      times=length(Files_IVF_fresh) )  
for(i in c(1:length(Sex_IVF_fresh)) ) {
  Sex_IVF_fresh[i] = "boy"
  if(i>=4) { Sex_IVF_fresh[i] = "boy" }
}



Files_IVF_frozen <- c(
  paste( inputDir,   "9_W1365C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "9_Q5-W1365D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "10_W1733C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "10_W1733D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "11_W1398C-boy-IVF-frozen_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir,   "11_W1398D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ) 
)
mySampleID_IVF_frozen <- rep( x="IVF_frozen", times=length(Files_IVF_frozen) )
for(i in c(1:length(mySampleID_IVF_frozen)) ) {
  mySampleID_IVF_frozen[i] = paste(mySampleID_IVF_frozen[i], i, sep="_")
}
myTreatment_IVF_frozen <- rep( x=2,  times=length(Files_IVF_frozen) )
Sex_IVF_frozen = rep( x="boy",      times=length(Files_IVF_frozen) )  
for(i in c(1:length(Sex_IVF_frozen)) ) {
  Sex_IVF_frozen[i] = "boy"
  if(i>=4) { Sex_IVF_frozen[i] = "boy" }
}



Files_ICSI_fresh <- c(
  paste( inputDir,   "12_W1579C-boy-ICSI-fresh-merge_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "12_Q17-W1579D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "13_W1647C-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "13_W1647D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "14_W1719C-boy-ICSI-fresh_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir,   "14_W1719D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
)
mySampleID_ICSI_fresh <- rep( x="ICSI_fresh", times=length(Files_ICSI_fresh) )
for(i in c(1:length(mySampleID_ICSI_fresh)) ) {
  mySampleID_ICSI_fresh[i] = paste(mySampleID_ICSI_fresh[i], i, sep="_")
}
myTreatment_ICSI_fresh <- rep( x=3,  times=length(Files_ICSI_fresh) )
Sex_ICSI_fresh = rep( x="boy",      times=length(Files_ICSI_fresh) )  
for(i in c(1:length(Sex_ICSI_fresh)) ) {
  Sex_ICSI_fresh[i] = "boy"
  if(i>=4) { Sex_ICSI_fresh[i] = "boy" }
}



Files_ICSI_frozen <- c(
  paste( inputDir,   "15_Q6-W871C-boy-ICSI-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir,   "15_Q4-W871D-boy-ICSI-frozen_Rep1.bismark.cov",    sep="/" ) 
)
mySampleID_ICSI_frozen <- rep( x="ICSI_frozen", times=length(Files_ICSI_frozen) )
for(i in c(1:length(mySampleID_ICSI_frozen)) ) {
  mySampleID_ICSI_frozen[i] = paste(mySampleID_ICSI_frozen[i], i, sep="_")
}
myTreatment_ICSI_frozen <- rep( x=4,  times=length(Files_ICSI_frozen) )
Sex_ICSI_frozen = rep( x="boy",      times=length(Files_ICSI_frozen) )  
for(i in c(1:length(Sex_ICSI_frozen)) ) {
  Sex_ICSI_frozen[i] = "boy"
  if(i>=4) { Sex_ICSI_frozen[i] = "boy" }
}




Files_All <- c(
  Files_NC,
  Files_IVF_fresh,
  Files_IVF_frozen,
  Files_ICSI_fresh,
  Files_ICSI_frozen  
)
Files_All2 <- as.list( Files_All )


mySampleID_All <- c( 
  mySampleID_NC,
  mySampleID_IVF_fresh,
  mySampleID_IVF_frozen,
  mySampleID_ICSI_fresh,
  mySampleID_ICSI_frozen  
)
mySampleID_All2 <- as.list( mySampleID_All )


myTreatment_All <- c( 
  myTreatment_NC,
  myTreatment_IVF_fresh,
  myTreatment_IVF_frozen,
  myTreatment_ICSI_fresh,
  myTreatment_ICSI_frozen 
)       
 




## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple

myType1_sex = c( 
  Sex_NC, Sex_IVF_fresh,   Sex_IVF_frozen,  Sex_ICSI_fresh,  Sex_ICSI_frozen
)

myType1_sex_shape  = myType1_sex
for(i in c(1:length(myType1_sex_shape)) ) {
  if(myType1_sex_shape[i] == "boy")    { myType1_sex_shape[i] = c("boy"=16) }
  if(myType1_sex_shape[i] == "girl")   { myType1_sex_shape[i] = c("girl"=17) }
  if(myType1_sex_shape[i] == "father") { myType1_sex_shape[i] = c("father"=1) }
  if(myType1_sex_shape[i] == "mother") { myType1_sex_shape[i] = c("mother"=2) }
}

myType1_sex_shape2 =  c( "boy"=16,  "girl"=17,  "father"=1, "mother"=2 )  




myType2_Tech = c( 
  rep( x="NC",               times=length(Files_NC) ) ,
  rep( x="IVF_fresh",        times=length(Files_IVF_fresh) ) ,
  rep( x="IVF_frozen",       times=length(Files_IVF_frozen) ) ,
  rep( x="ICSI_fresh",       times=length(Files_ICSI_fresh) ) ,
  rep( x="ICSI_frozen",      times=length(Files_ICSI_frozen) )  
)
myType2_Tech_color  = c( 
  rep( c("NC"="black"),                times=length(Files_NC) ) ,
  rep( c("IVF_fresh"="blue"),          times=length(Files_IVF_fresh) ) ,
  rep( c("IVF_frozen"="green"),        times=length(Files_IVF_frozen) ) ,
  rep( c("ICSI_fresh"="red"),          times=length(Files_ICSI_fresh) ) ,
  rep( c("ICSI_frozen"="purple"),      times=length(Files_ICSI_frozen) )  
)
myType2_Tech_color2 = c(
  "NC"="black", "IVF_fresh"="blue", "IVF_frozen"="green", 
  "ICSI_fresh"="red", "ICSI_frozen"="purple"
)
## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple





length( Files_All )
length( Files_All2 )
length( mySampleID_All )
length( mySampleID_All2 )
length( myTreatment_All )
length( myType1_sex )
length( myType1_sex_shape )
length( myType2_Tech )
length( myType2_Tech_color )






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




MyCluster_1 <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  dev.off()
}




MyCluster_2 <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="canberra",    method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  dev.off()
}







#################################################################
myOutDir_sub1 = paste(myOutDir, "/1-ReadRawFiles",  sep="");
if( ! file.exists(myOutDir_sub1) ) { dir.create(myOutDir_sub1, recursive = TRUE) }



sink( file=paste(myOutDir_sub1, "1-length-variables.txt", sep="/") )

length( Files_All )
length( Files_All2 )
length( mySampleID_All )
length( mySampleID_All2 )
length( myTreatment_All )
length( myType1_sex )
length( myType1_sex_shape )
length( myType2_Tech )
length( myType2_Tech_color )

cat("############################")
print( "myTreatment_All:" )
print( myTreatment_All )
cat("############################")
print( "mySampleID_All:" )
print( mySampleID_All )
cat("############################")
print( "mySampleID_All2:" )
print( mySampleID_All2 )
cat("############################")
print( "Files_All:" )
print( Files_All )
cat("############################")
print( "Files_All2:" )
print( Files_All2 )

sink()





# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1, "2-theLog-of-read-inputFiles.txt", sep="/") )
myobj = methRead( Files_All2,
               sample.id = mySampleID_All2,
               assembly  = "hg38",
               treatment = myTreatment_All,
               context   = "CpG",
               pipeline  = "bismarkCoverage",
               mincov    = 1,       ## >= n
               header    = FALSE
)
sink()




sink( file=paste(myOutDir_sub1, "3-all-rawFiles.txt", sep="/") )
    print(Files_All2)
    print("#########################")
    print("#########################")
    print(myobj)
sink()


 
sink( file=paste(myOutDir_sub1, "4-dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All2)) ) {
  print( "######################" )
  print(   Files_All2[[i]]  )
  print(   dim(myobj[[i]])  )
}
sink()



sink( file=paste(myOutDir_sub1, "5-dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All2)) ) {
  print(   dim(myobj[[i]])  )
}
sink()



continue_on_error <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error) 



myobj_nor <- normalizeCoverage(myobj)


######################################################################################################################################################


















MyCluster_3 <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyCluster_1(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    
  
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.03.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.03 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.08.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.08 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.12.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.12 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25)                    
  MyCluster_2(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30)                    
}







MyPCA_1 <- function(  mat2,   path2,   file2  ) {
  PCA1 <- prcomp( t(mat2)  )
  print( names(PCA1) )
  
  sink( file = paste(path2, "/1A_PCA_results_",  file2,  ".txt",  sep="") )
  print( PCA1 )
  print( names(PCA1) )
  sink()
  
  sink( file = paste(path2, "/1B_PCA_summary_",  file2,  ".txt",   sep="") )
  print( summary(PCA1) )
  sink()
  
  
  sink( file = paste(path2, "/2_PCA_all_",  file2,  ".txt",   sep="") )
  print("####################### PCA1$sdev #########################")
  print(PCA1$sdev)
  print("####################### PCA1$rotation #########################")
  print(PCA1$rotation)
  print("####################### PCA1$center #########################")
  print(PCA1$center)
  print("####################### PCA1$scale #########################")
  print(PCA1$scale)
  print("####################### PCA1$x #########################")
  print(PCA1$x)
  sink()
  
  
  pdf( file = paste(path2, "/3_PCA_info_",  file2,  ".pdf",   sep="")  )
  plot(PCA1, type="lines")
  print( fviz_eig(PCA1) )
  dev.off() 
  
  
  
  
  my_fviz_pca_ind1 <- fviz_pca_ind(PCA1,
                                                 col.ind = "cos2", # Color by the quality of representation
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE     # Avoid text overlapping
  )
  my_fviz_pca_ind2 <- fviz_pca_ind(PCA1,
                                                 col.ind =  as.factor( myTreatment_All ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 addEllipses = TRUE, # Concentration ellipses
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE
  )
  my_fviz_pca_ind3 <- fviz_pca_ind(PCA1,
                                                 col.ind =  as.factor( myTreatment_All ) , # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 addEllipses = TRUE, # Concentration ellipses
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE, 
                                                 label = "none", 
                                                 alpha.ind = 1
  )
  my_fviz_pca_ind4<- fviz_pca_ind(PCA1,
                                                 col.ind =  as.factor( myTreatment_All ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 #legend.title = "Groups",
                                                 repel = TRUE, 
                                                 label = "none", 
                                                 alpha.ind = 1
  )
  
  
  svg( file=paste(path2, "/4A-PCA-2D-1", file2, ".svg",  sep="") )
  print(my_fviz_pca_ind1)
  dev.off() 
  
  svg( file=paste(path2, "/4B-PCA-2D-2", file2, ".svg",  sep="")  )
  print(my_fviz_pca_ind2)
  dev.off() 
  
  svg( file=paste(path2, "/4C-PCA-2D-3", file2, ".svg",  sep="")  )
  print(my_fviz_pca_ind3)
  dev.off() 
  
  svg( file=paste(path2, "/4D-PCA-2D-4", file2, ".svg",  sep="")  )
  print(my_fviz_pca_ind4)
  dev.off() 
}








 

######################################################################################################################################################
######################################################################################################################################################
myOutDir_2two = paste(myOutDir, "/2A-NC-vs-IVFfresh-2kbBin",  sep="");
if( ! file.exists(myOutDir_2two) ) { dir.create(myOutDir_2two, recursive = TRUE) }



myobj_2two <- reorganize(myobj_nor, sample.ids = c(mySampleID_NC,  mySampleID_IVF_fresh), 
                      treatment  = c(myTreatment_NC, myTreatment_IVF_fresh)  
                      )
sink( file=paste(myOutDir_2two , "1-select-subSets.txt", sep="/")  )
length(myobj_2two) 
getSampleID(myobj_2two) 
getTreatment(myobj_2two) 
sink()


                  
tiles_2two = tileMethylCounts( myobj_2two,   win.size=2000,   step.size=2000) ## 2kb bin

meth_2two  = unite( tiles_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
mat_2two   = percMethylation( meth_2two )


sink( file=paste(myOutDir_2two , "0-dimensions-2kbBin.txt", sep="/")  )
print( tiles_2two_sub2_2kb )
print("#########dimensions:")
print( dim(meth_2two_sub2_2kb)  )   
print( dim(mat_2two_sub2_2kb)   )
sink()







write.table(meth_2two_sub2_2kb , 
            file = paste(myOutDir_2two,   "1A-meth-2kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two_sub2_2kb , 
            file = paste(myOutDir_2two,   "1B-mat-2kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





pdf( file=paste(myOutDir_2two, "2A-MethylationStats-2kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_2two_sub2_2kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_2two, "2B-MethylationStats-2kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_2two_sub2_2kb[[i]] )  )
}
sink()




pdf( file=paste(myOutDir_2two, "3A-CoverageStats-2kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_2two_sub2_2kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_2two, "3B-CoverageStats-2kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_2two_sub2_2kb[[i]] )  )
}
sink()
 

 



meth_2two_sub2_2kb_1 <- reorganize(meth_2two_sub2_2kb, 
                      sample.ids = c(mySampleID_NC, mySampleID_IVF_fresh), 
                      treatment  = c(myTreatment_NC, myTreatment_IVF_fresh)  
)
MyCluster_3(  mymeth2=meth_2two_sub2_2kb_1 ,  path2=myOutDir_2two,   
              file2="4-cluster", width2=8, height2=5 )




MyPCA_1(  mat2=mat_2two_sub2_2kb,   path2=myOutDir_2two,   file2="5-PCA"  )
  
  
PCASamples()



pdf( file=paste(myOutDir_2two, "5-PCA-2kbBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_2two_sub2_2kb , screeplot=TRUE)
PCASamples(meth_2two_sub2_2kb )
dev.off()

#####################








##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.

 
myres.PCA1  <- PCA(t(mat_2two_sub2_2kb),   graph = FALSE)

sink( file = paste(myOutDir_2two, "20_PCA-info-2kbBin-byPCA.txt",  sep="/") )
print( myres.PCA1 )
print( "#################################" )
print( summary(myres.PCA1) )
print( "#################################"  )
myeig.val_2two_sub2_2kb <- get_eigenvalue( myres.PCA1 )
myeig.val_2two_sub2_2kb
sink() 


pdf( file = paste(myOutDir_2two, "21_PCA-screePlot-2kbBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.PCA1, addlabels = TRUE )
fviz_screeplot(X=myres.PCA1, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.PCA1, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.PCA1, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.PCA1, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_2two_sub2_2kb <- fviz_pca_ind(myres.PCA1,
                                               col.ind = "cos2", # Color by the quality of representation
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_2two_sub2_2kb <- fviz_pca_ind(myres.PCA1,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)
my_fviz_pca_ind3_2two_sub2_2kb <- fviz_pca_ind(myres.PCA1,
                                               col.ind =  as.factor( myTreatment ) , # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind4_2two_sub2_2kb <- fviz_pca_ind(myres.PCA1,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               #legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind5_2two_sub2_2kb <- fviz_pca_ind(myres.PCA1,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)

svg(file=paste(myOutDir_2two, "22A-PCA-2D-1-2kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_2two_sub2_2kb)
dev.off() 

svg(file=paste(myOutDir_2two, "22B-PCA-2D-2-2kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_2two_sub2_2kb)
dev.off() 

svg(file=paste(myOutDir_2two, "22C-PCA-2D-3-2kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_2two_sub2_2kb)
dev.off() 

svg(file=paste(myOutDir_2two, "22D-PCA-2D-4-2kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_2two_sub2_2kb)
dev.off() 

svg(file=paste(myOutDir_2two, "22E-PCA-2D-4-2kbBin.svg", sep="/") )
print(my_fviz_pca_ind5_2two_sub2_2kb)
dev.off() 



#############################


myres.PCA1_matrix <- myres.PCA1$ind$coord
dim( myres.PCA1_matrix )

myres.PCA1_Contri  <- (myres.PCA1$eig)[,2] 
myres.PCA1_Contri  <- myres.PCA1_Contri/sum(myres.PCA1_Contri)
myres.PCA1_Contri  <- myres.PCA1_Contri * 100
myres.PCA1_Contri  <- round(myres.PCA1_Contri, 2)

label1_2two_sub2_2kb <-   paste( "PC1 ",  "(", myres.PCA1_Contri[1], "%)", sep="" )
label2_2two_sub2_2kb <-   paste( "PC2 ",  "(", myres.PCA1_Contri[2], "%)", sep="" )
label3_2two_sub2_2kb <-   paste( "PC3 ",  "(", myres.PCA1_Contri[3], "%)", sep="" )
label1_2two_sub2_2kb  
label2_2two_sub2_2kb  
label3_2two_sub2_2kb  


myLabel = myType1
dataframeB_2two_sub2_2kb  <- data.frame( as.data.frame(myres.PCA1_matrix), myType2, myType3, myLabel   ) 
dataframeB_2two_sub2_2kb 


FigureTemp1_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_2kb) +   ylab(label2_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="23A-PCA-PC1-PC2-2kbBin",  height1=3,  width1=4.8)


FigureTemp2_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_2kb) +   ylab(label2_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="23B-PCA-PC1-PC2-alpha-2kbBin",   height1=3,  width1=4.8)


FigureTemp3_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_2kb) +   ylab(label2_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="23C-PCA-PC1-PC2-smallDot-2kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_2kb) +   ylab(label2_2two_sub2_2kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="23D-PCA-PC1-PC2-big-2kbBin",   height1=3,  width1=4.8)



FigureTemp5_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_2kb) +   ylab(label2_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="23E-PCA-PC1-PC2-text-2kbBin",   height1=3,  width1=4.8)


FigureTemp6_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_2kb) +   ylab(label2_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="23F-PCA-PC1-PC2-text2-2kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_2two_sub2_2kb) +   ylab(label3_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="24A-PCA-PC1-PC3-2kbBin",   height1=3,  width1=4.8)


FigureTemp12_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_2two_sub2_2kb) +   ylab(label3_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="24B-PCA-PC1-PC3-alpha-2kbBin",   height1=3,  width1=4.8)


FigureTemp13_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_2two_sub2_2kb) +   ylab(label3_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="24C-PCA-PC1-PC3-smallDot-2kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_2two_sub2_2kb) +   ylab(label3_2two_sub2_2kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="24D-PCA-PC1-PC3-big-2kbBin",  height1=3,  width1=4.8)



FigureTemp15_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_2two_sub2_2kb) +   ylab(label3_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="24E-PCA-PC1-PC3-text-2kbBin",   height1=3,  width1=4.8)


FigureTemp16_2two_sub2_2kb  <- ggplot( data = dataframeB_2two_sub2_2kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_2two_sub2_2kb) +   ylab(label3_2two_sub2_2kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_2two_sub2_2kb ,  path1=myOutDir_2two, fileName1="24F-PCA-PC1-PC3-text2-2kbBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_2two, "25A_PCA-3d-2kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_2two_sub2_2kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_2kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_2two_sub2_2kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_2kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_2two_sub2_2kb[,1:3], pch =myType2_shape, 
  xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_2kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_2two, "25B_PCA-3d-2kbBin-withLabels.pdf",  sep="/") )

s3d1_2two_sub2_2kb <- scatterplot3d( 
  dataframeB_2two_sub2_2kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_2kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_2two_sub2_2kb$xyz.convert(dataframeB_2two_sub2_2kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_2two_sub2_2kb <- scatterplot3d( 
  dataframeB_2two_sub2_2kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_2kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_2two_sub2_2kb$xyz.convert(dataframeB_2two_sub2_2kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_2two_sub2_2kb <- scatterplot3d( 
  dataframeB_2two_sub2_2kb[,1:3], pch =myType2_shape, 
  xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_2two_sub2_2kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_2two_sub2_2kb$xyz.convert(dataframeB_2two_sub2_2kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_2two, "26A_PCA-3d-2kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_2two, "26B_PCA-3d-2kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_2two, "27_PCA-3d-2kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )



##
scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )



##
scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )


scatter3D(x = dataframeB_2two_sub2_2kb[,1], y = dataframeB_2two_sub2_2kb[,2], z = dataframeB_2two_sub2_2kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_2two_sub2_2kb,  ylab = label2_2two_sub2_2kb,  zlab = label3_2two_sub2_2kb  )

dev.off()






#################



km.res1_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=2, iter.max = 20 )
km.res2_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=3, iter.max = 20  )
km.res3_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=4, iter.max = 20  )
km.res4_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=5, iter.max = 20 )
km.res5_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=6, iter.max = 20 )
km.res6_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=7, iter.max = 20 )
km.res7_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=8, iter.max = 20 )
km.res8_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=9, iter.max = 20 )
km.res9_2two_sub2_2kb  <- kmeans( t(mat_2two_sub2_2kb), centers=10, iter.max = 20 )
km.res10_2two_sub2_2kb <- kmeans( t(mat_2two_sub2_2kb), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_2two, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_2two_sub2_2kb$cluster
print("#################### 3 classes:")
km.res2_2two_sub2_2kb$cluster
print("#################### 4 classes:")
km.res3_2two_sub2_2kb$cluster
print("#################### 5 classes:")
km.res4_2two_sub2_2kb$cluster
print("#################### 6 classes:")
km.res5_2two_sub2_2kb$cluster
print("#################### 7 classes:")
km.res6_2two_sub2_2kb$cluster
print("#################### 8 classes:")
km.res7_2two_sub2_2kb$cluster
print("#################### 9 classes:")
km.res8_2two_sub2_2kb$cluster
print("#################### 10 classes:")
km.res9_2two_sub2_2kb$cluster
print("#################### 11 classes:")
km.res10_2two_sub2_2kb$cluster
sink()



pdf( file = paste(myOutDir_2two, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_2kb), col = km.res1_2two_sub2_2kb$cluster )

plot( t(mat_2two_sub2_2kb), col = km.res1_2two_sub2_2kb$cluster)
points(km.res1_2two_sub2_2kb$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_2kb), col = km.res1_2two_sub2_2kb$cluster)
points(km.res1_2two_sub2_2kb$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_2two_sub2_2kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 
 

 



pdf( file = paste(myOutDir_2two, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_2kb), col = km.res2_2two_sub2_2kb$cluster )

plot( t(mat_2two_sub2_2kb), col = km.res2_2two_sub2_2kb$cluster)
points(km.res2_2two_sub2_2kb$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_2kb), col = km.res2_2two_sub2_2kb$cluster)
points(km.res2_2two_sub2_2kb$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_2two_sub2_2kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_2two, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_2kb), col = km.res3_2two_sub2_2kb$cluster )

plot( t(mat_2two_sub2_2kb), col = km.res3_2two_sub2_2kb$cluster)
points(km.res3_2two_sub2_2kb$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_2kb), col = km.res3_2two_sub2_2kb$cluster)
points(km.res3_2two_sub2_2kb$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_2two_sub2_2kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 


 




pdf( file = paste(myOutDir_2two, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_2two_sub2_2kb), col = km.res4_2two_sub2_2kb$cluster )

plot( t(mat_2two_sub2_2kb), col = km.res4_2two_sub2_2kb$cluster)
points(km.res4_2two_sub2_2kb$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_2two_sub2_2kb), col = km.res4_2two_sub2_2kb$cluster)
points(km.res4_2two_sub2_2kb$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_2two_sub2_2kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "euclidean")
res.dist2_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "maximum")
res.dist3_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "manhattan")
res.dist4_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "canberra")
res.dist5_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "binary")
res.dist6_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "minkowski")
res.dist7_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "pearson")
res.dist8_2two_sub2_2kb <- get_dist( t(mat_2two_sub2_2kb) ,   method = "spearman")
 



pdf( file = paste(myOutDir_2two, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_2two_sub2_2kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_2two, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_2two_sub2_2kb ) 
print("################### maximum: ")
print(res.dist2_2two_sub2_2kb ) 
print("################### manhattan: ")
print(res.dist3_2two_sub2_2kb ) 
print("################### canberra: ")
print(res.dist4_2two_sub2_2kb ) 
print("################### binary: ")
print(res.dist5_2two_sub2_2kb ) 
print("################### minkowski: ")
print(res.dist6_2two_sub2_2kb ) 
print("################### pearson: ")
print(res.dist7_2two_sub2_2kb ) 
print("################### spearman: ")
print(res.dist8_2two_sub2_2kb ) 
sink() 


# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "ward.D" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "ward.D" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "ward.D" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "ward.D" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "ward.D" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "ward.D" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "ward.D" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "ward.D" )  

pdf( file = paste(myOutDir_2two, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "ward.D2" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "ward.D2" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "ward.D2" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "ward.D2" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "ward.D2" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "ward.D2" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "ward.D2" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "ward.D2" )  

pdf( file = paste(myOutDir_2two, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "single" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "single" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "single" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "single" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "single" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "single" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "single" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "single" )  

pdf( file = paste(myOutDir_2two, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "complete" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "complete" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "complete" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "complete" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "complete" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "complete" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "complete" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "complete" )  

pdf( file = paste(myOutDir_2two, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "average" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "average" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "average" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "average" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "average" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "average" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "average" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "average" )  

pdf( file = paste(myOutDir_2two, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "mcquitty" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "mcquitty" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "mcquitty" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "mcquitty" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "mcquitty" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "mcquitty" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "mcquitty" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "mcquitty" )  

pdf( file = paste(myOutDir_2two, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "median" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "median" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "median" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "median" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "median" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "median" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "median" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "median" )  

pdf( file = paste(myOutDir_2two, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_2two_sub2_2kb <- hclust(res.dist1_2two_sub2_2kb, method = "centroid" )   
res.hc2_2two_sub2_2kb <- hclust(res.dist2_2two_sub2_2kb, method = "centroid" )   
res.hc3_2two_sub2_2kb <- hclust(res.dist3_2two_sub2_2kb, method = "centroid" )   
res.hc4_2two_sub2_2kb <- hclust(res.dist4_2two_sub2_2kb, method = "centroid" )   
res.hc5_2two_sub2_2kb <- hclust(res.dist5_2two_sub2_2kb, method = "centroid" )   
res.hc6_2two_sub2_2kb <- hclust(res.dist6_2two_sub2_2kb, method = "centroid" )   
res.hc7_2two_sub2_2kb <- hclust(res.dist7_2two_sub2_2kb, method = "centroid" )   
res.hc8_2two_sub2_2kb <- hclust(res.dist8_2two_sub2_2kb, method = "centroid" )  

pdf( file = paste(myOutDir_2two, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_2two_sub2_2kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_2two, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb )
fviz_dend(res.hc2_2two_sub2_2kb )
fviz_dend(res.hc3_2two_sub2_2kb )
fviz_dend(res.hc4_2two_sub2_2kb )
fviz_dend(res.hc5_2two_sub2_2kb )
fviz_dend(res.hc6_2two_sub2_2kb )
fviz_dend(res.hc7_2two_sub2_2kb )
fviz_dend(res.hc8_2two_sub2_2kb )
dev.off() 


pdf( file = paste(myOutDir_2two, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc2_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc3_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc4_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc5_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc6_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc7_2two_sub2_2kb, label_cols = myType3_color )
fviz_dend(res.hc8_2two_sub2_2kb, label_cols = myType3_color )
dev.off() 









##################
myOutDir_2two_diffMe   <- paste(myOutDir_2two, "diffMe_DMC_DMR",  sep="/")
if( ! file.exists(myOutDir_2two_diffMe) ) { dir.create(myOutDir_2two_diffMe, recursive = TRUE) }


myDiff_2two_sub2_2kb = calculateDiffMeth(meth_2two_sub2_2kb,  num.cores=16)
dim(myDiff_2two_sub2_2kb)
names(myDiff_2two_sub2_2kb)
head(myDiff_2two_sub2_2kb)

myQvalue_2two_sub2_2kb = myDiff_2two_sub2_2kb$qvalue
myMethDi_2two_sub2_2kb = myDiff_2two_sub2_2kb$meth.diff



pdf(paste(myOutDir_2two_diffMe,  "1A_qvalue_distribution.pdf",  sep="/"))
hist( myQvalue_2two_sub2_2kb,  nclass=100, xlim=c(0, 1), freq=FALSE)
hist( myQvalue_2two_sub2_2kb[myQvalue_2two_sub2_2kb<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE)
hist( myQvalue_2two_sub2_2kb,  nclass=100, xlim=c(0, 1), freq=TRUE)
hist( myQvalue_2two_sub2_2kb[myQvalue_2two_sub2_2kb<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE)
dev.off() 


pdf(paste(myOutDir_2two_diffMe,  "1B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_2two_sub2_2kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_2two_sub2_2kb, nclass=100, xlim=c(0, 50), freq=FALSE)
hist(myMethDi_2two_sub2_2kb[myMethDi_2two_sub2_2kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_2two_sub2_2kb[myMethDi_2two_sub2_2kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_2two_sub2_2kb[myMethDi_2two_sub2_2kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
hist(myMethDi_2two_sub2_2kb, nclass=100, xlim=c(0, 100), freq=TRUE)
hist(myMethDi_2two_sub2_2kb, nclass=100, xlim=c(0, 50), freq=TRUE)
hist(myMethDi_2two_sub2_2kb[myMethDi_2two_sub2_2kb>=20], nclass=81, xlim=c(20, 100),  freq=TRUE)
hist(myMethDi_2two_sub2_2kb[myMethDi_2two_sub2_2kb>=10], nclass=91, xlim=c(10, 100),  freq=TRUE)
hist(myMethDi_2two_sub2_2kb[myMethDi_2two_sub2_2kb<=10], nclass=100, xlim=c(0, 10),   freq=TRUE)
dev.off() 




myColor1_2two_sub2_2kb <- rep( "no",   times= length(myQvalue_2two_sub2_2kb) )
myColor1_2two_sub2_2kb[ (abs(myMethDi_2two_sub2_2kb)>10) & (myQvalue_2two_sub2_2kb<0.001) ]  <- "yes"
length( myColor1_2two_sub2_2kb[myColor1_2two_sub2_2kb=="yes"] )



DataFrame2_2two_sub2_2kb <- data.frame(myx1 =  myMethDi_2two_sub2_2kb,   
                                       myy1 =  -log10(myQvalue_2two_sub2_2kb),  
                                       mycolor1 = myColor1_2two_sub2_2kb )

ggplot(DataFrame2_2two_sub2_2kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_diffMe,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_2two_sub2_2kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_2two_diffMe,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_2two_sub2_2kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_2two_diffMe,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )




ggplot(DataFrame2_2two_sub2_2kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_2two_diffMe,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_2two_sub2_2kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_2two_diffMe,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_2two_sub2_2kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_2two_diffMe,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )









myDiff25p.hypo_2two_sub2_2kb  = getMethylDiff(myDiff_2two_sub2_2kb, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_2two_sub2_2kb = getMethylDiff(myDiff_2two_sub2_2kb, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
myDiff25p_2two_sub2_2kb       = getMethylDiff(myDiff_2two_sub2_2kb, difference=10, qvalue=0.05)
myDiffTemp_2two_sub2_2kb      = getMethylDiff(myDiff_2two_sub2_2kb, difference=0,  qvalue=0.05)

write.table(myDiff_2two_sub2_2kb , file = paste(myOutDir_2two_diffMe,"2A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_2two_sub2_2kb , file = paste(myOutDir_2two_diffMe,"2B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_2two_sub2_2kb , file = paste(myOutDir_2two_diffMe,"2C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_2two_sub2_2kb , file = paste(myOutDir_2two_diffMe,"2D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_2two_sub2_2kb , file = paste(myOutDir_2two_diffMe,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_2two_diffMe, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_2two_sub2_2kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_2two_diffMe, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_2two_sub2_2kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
diffGeneAnn_2two_sub2_2kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub2_2kb,  "GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_diffMe, "4A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffGeneAnn_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub2_2kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2_2kb,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "4C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffCpGann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2_2kb,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "4E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffrepeatann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_2kb, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "5A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_2kb, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "5C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_2kb, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "5E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_2kb, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "5G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub2_2kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2_2kb, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "5I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_2kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub2_2kb_hypo)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub2_2kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_2two_sub2_2kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub2_2kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_2two_diffMe, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffGeneAnn_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_2two_sub2_2kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2_2kb,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffCpGann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2_2kb,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffrepeatann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_2kb, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_2kb, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_2kb, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_2kb, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_2two_sub2_2kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2_2kb, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_2two_diffMe, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_2kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_2two_sub2_2kb_hyper)
sink()

pdf( file=paste(myOutDir_2two_diffMe, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_2two_sub2_2kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##################################################################################################################
##################################################################################################################







































######################################################################################################################################################
######################################################################################################################################################
myOutDir_3three = paste(myOutDir, "/3-1kbBin",  sep="");
if( ! file.exists(myOutDir_3three) ) { dir.create(myOutDir_3three, recursive = TRUE) }


#####
tiles_3three_sub3_1kb = tileMethylCounts( myobj, win.size=1000, step.size=1000) ## 1kb bin
meth_3three_sub3_1kb  = unite( tiles_3three_sub3_1kb, destrand=FALSE   )   ## 100% overlap
mat_3three_sub3_1kb   = percMethylation( meth_3three_sub3_1kb )


sink( file=paste(myOutDir_3three , "0-dimensions-1kbBin.txt", sep="/")  )
print( tiles_3three_sub3_1kb )
print("#########dimensions:")
print( dim(meth_3three_sub3_1kb)  )   
print( dim(mat_3three_sub3_1kb)   )
sink()




write.table(meth_3three_sub3_1kb , 
            file = paste(myOutDir_3three,   "1A-meth-1kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_3three_sub3_1kb , 
            file = paste(myOutDir_3three,   "1B-mat-1kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





pdf( file=paste(myOutDir_3three, "2A-MethylationStats-1kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_3three_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_3three, "2B-MethylationStats-1kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_3three_sub3_1kb[[i]] )  )
}
sink()




pdf( file=paste(myOutDir_3three, "3A-CoverageStats-1kbBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_3three_sub3_1kb[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_3three, "3B-CoverageStats-1kbBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_3three_sub3_1kb[[i]] )  )
}
sink()





MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4A-cluster-kept100percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0 )

MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4B-cluster-kept80percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.2 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4C-cluster-kept60percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.4 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4D-cluster-kept50percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.5 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4E-cluster-kept40percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.6 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4F-cluster-kept30percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.7 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4G-cluster-kept20percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.8 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4H-cluster-kept10percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.9 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4I-cluster-kept5percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.95 )


MyCluster_1(  mymeth1 = meth_3three_sub3_1kb , path1 = myOutDir_3three, 
              file1 = "4J-cluster-kept1percent-1kbBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.99 )





pdf( file=paste(myOutDir_3three, "5-PCA-1kbBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_3three_sub3_1kb , screeplot=TRUE)
PCASamples(meth_3three_sub3_1kb )
dev.off()

#####################





PCA_3three_sub3_1kb <- prcomp( t(mat_3three_sub3_1kb)  )
names(PCA_3three_sub3_1kb)

sink( file = paste(myOutDir_3three,  "6A-PCA-1kbBin.txt",  sep="/") )
print(PCA_3three_sub3_1kb)
sink()

sink( file = paste(myOutDir_3three,  "6B-PCA-summary-1kbBin.txt",  sep="/") )
summary(PCA_3three_sub3_1kb)
sink()


sink( file = paste(myOutDir_3three,  "6C-PCA-all-1kbBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_3three_sub3_1kb$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_3three_sub3_1kb$rotation)
print("####################### myPCA2$center #########################")
print(PCA_3three_sub3_1kb$center)
print("####################### myPCA2$scale #########################")
print(PCA_3three_sub3_1kb$scale)
print("####################### myPCA2$x #########################")
print(PCA_3three_sub3_1kb$x)
sink()


pdf( file=paste(myOutDir_3three,   "6D-PCA-info-1kbBin.pdf", sep="/")  )
plot(PCA_3three_sub3_1kb, type="lines")
fviz_eig(PCA_3three_sub3_1kb)
dev.off() 




my_fviz_pca_ind1_3three_sub3_1kb <- fviz_pca_ind(PCA_3three_sub3_1kb,
                                                 col.ind = "cos2", # Color by the quality of representation
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_3three_sub3_1kb <- fviz_pca_ind(PCA_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 addEllipses = TRUE, # Concentration ellipses
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE
)
my_fviz_pca_ind3_3three_sub3_1kb <- fviz_pca_ind(PCA_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ) , # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 addEllipses = TRUE, # Concentration ellipses
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE, 
                                                 label = "none", 
                                                 alpha.ind = 1
)
my_fviz_pca_ind4_3three_sub3_1kb <- fviz_pca_ind(PCA_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 #legend.title = "Groups",
                                                 repel = TRUE, 
                                                 label = "none", 
                                                 alpha.ind = 1
)


svg(file=paste(myOutDir_3three, "7A-PCA-2D-1-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "7B-PCA-2D-2-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "7C-PCA-2D-3-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "7D-PCA-2D-4-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_3three_sub3_1kb)
dev.off() 








#############################


PCA_3three_sub3_1kb_matrix <- PCA_3three_sub3_1kb$x
dim( PCA_3three_sub3_1kb_matrix )

PCA_3three_sub3_1kb_Contri  <- (PCA_3three_sub3_1kb$sdev)^2
PCA_3three_sub3_1kb_Contri  <- PCA_3three_sub3_1kb_Contri/sum(PCA_3three_sub3_1kb_Contri)
PCA_3three_sub3_1kb_Contri  <- PCA_3three_sub3_1kb_Contri * 100
PCA_3three_sub3_1kb_Contri  <- round(PCA_3three_sub3_1kb_Contri, 2)

label1_3three_sub3_1kb <-   paste( "PC1 ",  "(", PCA_3three_sub3_1kb_Contri[1], "%)", sep="" )
label2_3three_sub3_1kb <-   paste( "PC2 ",  "(", PCA_3three_sub3_1kb_Contri[2], "%)", sep="" )
label3_3three_sub3_1kb <-   paste( "PC3 ",  "(", PCA_3three_sub3_1kb_Contri[3], "%)", sep="" )
label1_3three_sub3_1kb  
label2_3three_sub3_1kb  
label3_3three_sub3_1kb  


myLabel = myType1
dataframeA_3three_sub3_1kb  <- data.frame( as.data.frame(PCA_3three_sub3_1kb_matrix), myType2, myType3, myLabel   ) 
dataframeA_3three_sub3_1kb 


FigureTemp1_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="8A-PCA-PC1-PC2-1kbBin",  height1=3,  width1=4.8)


FigureTemp2_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="8B-PCA-PC1-PC2-alpha-1kbBin",   height1=3,  width1=4.8)


FigureTemp3_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="8C-PCA-PC1-PC2-smallDot-1kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="8D-PCA-PC1-PC2-big-1kbBin",   height1=3,  width1=4.8)



FigureTemp5_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="8E-PCA-PC1-PC2-text-1kbBin",   height1=3,  width1=4.8)


FigureTemp6_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="8F-PCA-PC1-PC2-text2-1kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="9A-PCA-PC1-PC3-1kbBin",   height1=3,  width1=4.8)


FigureTemp12_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="9B-PCA-PC1-PC3-alpha-1kbBin",   height1=3,  width1=4.8)


FigureTemp13_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="9C-PCA-PC1-PC3-smallDot-1kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="9D-PCA-PC1-PC3-big-1kbBin",  height1=3,  width1=4.8)



FigureTemp15_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="9E-PCA-PC1-PC3-text-1kbBin",   height1=3,  width1=4.8)


FigureTemp16_3three_sub3_1kb  <- ggplot( data = dataframeA_3three_sub3_1kb, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="9F-PCA-PC1-PC3-text2-1kbBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_3three, "10A_PCA-3d-1kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_3three_sub3_1kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_1kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_3three_sub3_1kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_1kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_3three_sub3_1kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_1kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_3three, "10B_PCA-3d-1kbBin-withLabels.pdf",  sep="/") )

s3d1_3three_sub3_1kb <- scatterplot3d( 
  dataframeA_3three_sub3_1kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_1kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_3three_sub3_1kb$xyz.convert(dataframeA_3three_sub3_1kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_3three_sub3_1kb <- scatterplot3d( 
  dataframeA_3three_sub3_1kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_1kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_3three_sub3_1kb$xyz.convert(dataframeA_3three_sub3_1kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_3three_sub3_1kb <- scatterplot3d( 
  dataframeA_3three_sub3_1kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_3three_sub3_1kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_3three_sub3_1kb$xyz.convert(dataframeA_3three_sub3_1kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "11A_PCA-3d-1kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "11B_PCA-3d-1kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_3three, "12_PCA-3d-1kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )



##
scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )



##
scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeA_3three_sub3_1kb[,1], y = dataframeA_3three_sub3_1kb[,2], z = dataframeA_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.


myres.pca_3three_sub3_1kb  <- PCA(t(mat_3three_sub3_1kb),   graph = FALSE)

sink( file = paste(myOutDir_3three, "20_PCA-info-1kbBin-byPCA.txt",  sep="/") )
print( myres.pca_3three_sub3_1kb )
print( "#################################" )
print( summary(myres.pca_3three_sub3_1kb) )
print( "#################################"  )
myeig.val_3three_sub3_1kb <- get_eigenvalue( myres.pca_3three_sub3_1kb )
myeig.val_3three_sub3_1kb
sink() 


pdf( file = paste(myOutDir_3three, "21_PCA-screePlot-1kbBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_3three_sub3_1kb, addlabels = TRUE )
fviz_screeplot(X=myres.pca_3three_sub3_1kb, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_3three_sub3_1kb, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_3three_sub3_1kb, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_3three_sub3_1kb, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_3three_sub3_1kb <- fviz_pca_ind(myres.pca_3three_sub3_1kb,
                                                 col.ind = "cos2", # Color by the quality of representation
                                                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                 repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_3three_sub3_1kb <- fviz_pca_ind(myres.pca_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 addEllipses = TRUE, # Concentration ellipses
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE
)
my_fviz_pca_ind3_3three_sub3_1kb <- fviz_pca_ind(myres.pca_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ) , # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 addEllipses = TRUE, # Concentration ellipses
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE, 
                                                 label = "none", 
                                                 alpha.ind = 1
)
my_fviz_pca_ind4_3three_sub3_1kb <- fviz_pca_ind(myres.pca_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 #legend.title = "Groups",
                                                 repel = TRUE, 
                                                 label = "none", 
                                                 alpha.ind = 1
)
my_fviz_pca_ind5_3three_sub3_1kb <- fviz_pca_ind(myres.pca_3three_sub3_1kb,
                                                 col.ind =  as.factor( myTreatment ), # color by groups
                                                 #palette = c("#00AFBB",  "#FC4E07"),
                                                 ellipse.type = "confidence",
                                                 legend.title = "Groups",
                                                 repel = TRUE
)

svg(file=paste(myOutDir_3three, "22A-PCA-2D-1-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind1_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22B-PCA-2D-2-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind2_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22C-PCA-2D-3-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind3_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22D-PCA-2D-4-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind4_3three_sub3_1kb)
dev.off() 

svg(file=paste(myOutDir_3three, "22E-PCA-2D-4-1kbBin.svg", sep="/") )
print(my_fviz_pca_ind5_3three_sub3_1kb)
dev.off() 



#############################


myres.pca_3three_sub3_1kb_matrix <- myres.pca_3three_sub3_1kb$ind$coord
dim( myres.pca_3three_sub3_1kb_matrix )

myres.pca_3three_sub3_1kb_Contri  <- (myres.pca_3three_sub3_1kb$eig)[,2] 
myres.pca_3three_sub3_1kb_Contri  <- myres.pca_3three_sub3_1kb_Contri/sum(myres.pca_3three_sub3_1kb_Contri)
myres.pca_3three_sub3_1kb_Contri  <- myres.pca_3three_sub3_1kb_Contri * 100
myres.pca_3three_sub3_1kb_Contri  <- round(myres.pca_3three_sub3_1kb_Contri, 2)

label1_3three_sub3_1kb <-   paste( "PC1 ",  "(", myres.pca_3three_sub3_1kb_Contri[1], "%)", sep="" )
label2_3three_sub3_1kb <-   paste( "PC2 ",  "(", myres.pca_3three_sub3_1kb_Contri[2], "%)", sep="" )
label3_3three_sub3_1kb <-   paste( "PC3 ",  "(", myres.pca_3three_sub3_1kb_Contri[3], "%)", sep="" )
label1_3three_sub3_1kb  
label2_3three_sub3_1kb  
label3_3three_sub3_1kb  


myLabel = myType1
dataframeB_3three_sub3_1kb  <- data.frame( as.data.frame(myres.pca_3three_sub3_1kb_matrix), myType2, myType3, myLabel   ) 
dataframeB_3three_sub3_1kb 


FigureTemp1_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="23A-PCA-PC1-PC2-1kbBin",  height1=3,  width1=4.8)


FigureTemp2_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="23B-PCA-PC1-PC2-alpha-1kbBin",   height1=3,  width1=4.8)


FigureTemp3_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="23C-PCA-PC1-PC2-smallDot-1kbBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="23D-PCA-PC1-PC2-big-1kbBin",   height1=3,  width1=4.8)



FigureTemp5_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="23E-PCA-PC1-PC2-text-1kbBin",   height1=3,  width1=4.8)


FigureTemp6_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label2_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="23F-PCA-PC1-PC2-text2-1kbBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="24A-PCA-PC1-PC3-1kbBin",   height1=3,  width1=4.8)


FigureTemp12_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="24B-PCA-PC1-PC3-alpha-1kbBin",   height1=3,  width1=4.8)


FigureTemp13_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="24C-PCA-PC1-PC3-smallDot-1kbBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="24D-PCA-PC1-PC3-big-1kbBin",  height1=3,  width1=4.8)



FigureTemp15_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="24E-PCA-PC1-PC3-text-1kbBin",   height1=3,  width1=4.8)


FigureTemp16_3three_sub3_1kb  <- ggplot( data = dataframeB_3three_sub3_1kb, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_3three_sub3_1kb) +   ylab(label3_3three_sub3_1kb) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_3three_sub3_1kb ,  path1=myOutDir_3three, fileName1="24F-PCA-PC1-PC3-text2-1kbBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_3three, "25A_PCA-3d-1kbBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_3three_sub3_1kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_1kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_3three_sub3_1kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_1kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_3three_sub3_1kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_1kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_3three, "25B_PCA-3d-1kbBin-withLabels.pdf",  sep="/") )

s3d1_3three_sub3_1kb <- scatterplot3d( 
  dataframeB_3three_sub3_1kb[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_1kb$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_3three_sub3_1kb$xyz.convert(dataframeB_3three_sub3_1kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_3three_sub3_1kb <- scatterplot3d( 
  dataframeB_3three_sub3_1kb[,1:3], color=myType3_color, pch=19,   
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_1kb$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_3three_sub3_1kb$xyz.convert(dataframeB_3three_sub3_1kb[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_3three_sub3_1kb <- scatterplot3d( 
  dataframeB_3three_sub3_1kb[,1:3], pch =myType2_shape, 
  xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_3three_sub3_1kb$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_3three_sub3_1kb$xyz.convert(dataframeB_3three_sub3_1kb[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "26A_PCA-3d-1kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_3three, "26B_PCA-3d-1kbBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_3three, "27_PCA-3d-1kbBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )



##
scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )



##
scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )


scatter3D(x = dataframeB_3three_sub3_1kb[,1], y = dataframeB_3three_sub3_1kb[,2], z = dataframeB_3three_sub3_1kb[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_3three_sub3_1kb,  ylab = label2_3three_sub3_1kb,  zlab = label3_3three_sub3_1kb  )

dev.off()






#################



km.res1_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=2, iter.max = 20 )
km.res2_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=3, iter.max = 20  )
km.res3_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=4, iter.max = 20  )
km.res4_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=5, iter.max = 20 )
km.res5_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=6, iter.max = 20 )
km.res6_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=7, iter.max = 20 )
km.res7_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=8, iter.max = 20 )
km.res8_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=9, iter.max = 20 )
km.res9_3three_sub3_1kb  <- kmeans( t(mat_3three_sub3_1kb), centers=10, iter.max = 20 )
km.res10_3three_sub3_1kb <- kmeans( t(mat_3three_sub3_1kb), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_3three, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_3three_sub3_1kb$cluster
print("#################### 3 classes:")
km.res2_3three_sub3_1kb$cluster
print("#################### 4 classes:")
km.res3_3three_sub3_1kb$cluster
print("#################### 5 classes:")
km.res4_3three_sub3_1kb$cluster
print("#################### 6 classes:")
km.res5_3three_sub3_1kb$cluster
print("#################### 7 classes:")
km.res6_3three_sub3_1kb$cluster
print("#################### 8 classes:")
km.res7_3three_sub3_1kb$cluster
print("#################### 9 classes:")
km.res8_3three_sub3_1kb$cluster
print("#################### 10 classes:")
km.res9_3three_sub3_1kb$cluster
print("#################### 11 classes:")
km.res10_3three_sub3_1kb$cluster
sink()



pdf( file = paste(myOutDir_3three, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_1kb), col = km.res1_3three_sub3_1kb$cluster )

plot( t(mat_3three_sub3_1kb), col = km.res1_3three_sub3_1kb$cluster)
points(km.res1_3three_sub3_1kb$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_1kb), col = km.res1_3three_sub3_1kb$cluster)
points(km.res1_3three_sub3_1kb$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_3three_sub3_1kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 






pdf( file = paste(myOutDir_3three, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_1kb), col = km.res2_3three_sub3_1kb$cluster )

plot( t(mat_3three_sub3_1kb), col = km.res2_3three_sub3_1kb$cluster)
points(km.res2_3three_sub3_1kb$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_1kb), col = km.res2_3three_sub3_1kb$cluster)
points(km.res2_3three_sub3_1kb$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_3three_sub3_1kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_3three, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_1kb), col = km.res3_3three_sub3_1kb$cluster )

plot( t(mat_3three_sub3_1kb), col = km.res3_3three_sub3_1kb$cluster)
points(km.res3_3three_sub3_1kb$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_1kb), col = km.res3_3three_sub3_1kb$cluster)
points(km.res3_3three_sub3_1kb$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_3three_sub3_1kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 







pdf( file = paste(myOutDir_3three, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_3three_sub3_1kb), col = km.res4_3three_sub3_1kb$cluster )

plot( t(mat_3three_sub3_1kb), col = km.res4_3three_sub3_1kb$cluster)
points(km.res4_3three_sub3_1kb$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_3three_sub3_1kb), col = km.res4_3three_sub3_1kb$cluster)
points(km.res4_3three_sub3_1kb$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_3three_sub3_1kb), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "euclidean")
res.dist2_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "maximum")
res.dist3_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "manhattan")
res.dist4_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "canberra")
res.dist5_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "binary")
res.dist6_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "minkowski")
res.dist7_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "pearson")
res.dist8_3three_sub3_1kb <- get_dist( t(mat_3three_sub3_1kb) ,   method = "spearman")




pdf( file = paste(myOutDir_3three, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_3three_sub3_1kb,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_3three, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_3three_sub3_1kb ) 
print("################### maximum: ")
print(res.dist2_3three_sub3_1kb ) 
print("################### manhattan: ")
print(res.dist3_3three_sub3_1kb ) 
print("################### canberra: ")
print(res.dist4_3three_sub3_1kb ) 
print("################### binary: ")
print(res.dist5_3three_sub3_1kb ) 
print("################### minkowski: ")
print(res.dist6_3three_sub3_1kb ) 
print("################### pearson: ")
print(res.dist7_3three_sub3_1kb ) 
print("################### spearman: ")
print(res.dist8_3three_sub3_1kb ) 
sink() 


# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "ward.D" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "ward.D" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "ward.D" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "ward.D" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "ward.D" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "ward.D" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "ward.D" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "ward.D" )  

pdf( file = paste(myOutDir_3three, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "ward.D2" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "ward.D2" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "ward.D2" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "ward.D2" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "ward.D2" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "ward.D2" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "ward.D2" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "ward.D2" )  

pdf( file = paste(myOutDir_3three, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "single" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "single" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "single" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "single" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "single" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "single" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "single" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "single" )  

pdf( file = paste(myOutDir_3three, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "complete" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "complete" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "complete" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "complete" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "complete" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "complete" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "complete" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "complete" )  

pdf( file = paste(myOutDir_3three, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "average" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "average" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "average" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "average" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "average" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "average" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "average" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "average" )  

pdf( file = paste(myOutDir_3three, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "mcquitty" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "mcquitty" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "mcquitty" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "mcquitty" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "mcquitty" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "mcquitty" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "mcquitty" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "mcquitty" )  

pdf( file = paste(myOutDir_3three, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "median" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "median" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "median" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "median" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "median" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "median" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "median" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "median" )  

pdf( file = paste(myOutDir_3three, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_3three_sub3_1kb <- hclust(res.dist1_3three_sub3_1kb, method = "centroid" )   
res.hc2_3three_sub3_1kb <- hclust(res.dist2_3three_sub3_1kb, method = "centroid" )   
res.hc3_3three_sub3_1kb <- hclust(res.dist3_3three_sub3_1kb, method = "centroid" )   
res.hc4_3three_sub3_1kb <- hclust(res.dist4_3three_sub3_1kb, method = "centroid" )   
res.hc5_3three_sub3_1kb <- hclust(res.dist5_3three_sub3_1kb, method = "centroid" )   
res.hc6_3three_sub3_1kb <- hclust(res.dist6_3three_sub3_1kb, method = "centroid" )   
res.hc7_3three_sub3_1kb <- hclust(res.dist7_3three_sub3_1kb, method = "centroid" )   
res.hc8_3three_sub3_1kb <- hclust(res.dist8_3three_sub3_1kb, method = "centroid" )  

pdf( file = paste(myOutDir_3three, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_3three_sub3_1kb, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_3three, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb )
fviz_dend(res.hc2_3three_sub3_1kb )
fviz_dend(res.hc3_3three_sub3_1kb )
fviz_dend(res.hc4_3three_sub3_1kb )
fviz_dend(res.hc5_3three_sub3_1kb )
fviz_dend(res.hc6_3three_sub3_1kb )
fviz_dend(res.hc7_3three_sub3_1kb )
fviz_dend(res.hc8_3three_sub3_1kb )
dev.off() 


pdf( file = paste(myOutDir_3three, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc2_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc3_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc4_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc5_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc6_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc7_3three_sub3_1kb, label_cols = myType3_color )
fviz_dend(res.hc8_3three_sub3_1kb, label_cols = myType3_color )
dev.off() 









##################
myOutDir_3three_diffMe   <- paste(myOutDir_3three, "diffMe_DMC_DMR",  sep="/")
if( ! file.exists(myOutDir_3three_diffMe) ) { dir.create(myOutDir_3three_diffMe, recursive = TRUE) }


myDiff_3three_sub3_1kb = calculateDiffMeth(meth_3three_sub3_1kb,  num.cores=16)
dim(myDiff_3three_sub3_1kb)
names(myDiff_3three_sub3_1kb)
head(myDiff_3three_sub3_1kb)

myQvalue_3three_sub3_1kb = myDiff_3three_sub3_1kb$qvalue
myMethDi_3three_sub3_1kb = myDiff_3three_sub3_1kb$meth.diff



pdf(paste(myOutDir_3three_diffMe,  "1A_qvalue_distribution.pdf",  sep="/"))
hist( myQvalue_3three_sub3_1kb,  nclass=100, xlim=c(0, 1), freq=FALSE)
hist( myQvalue_3three_sub3_1kb[myQvalue_3three_sub3_1kb<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE)
hist( myQvalue_3three_sub3_1kb,  nclass=100, xlim=c(0, 1), freq=TRUE)
hist( myQvalue_3three_sub3_1kb[myQvalue_3three_sub3_1kb<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE)
dev.off() 


pdf(paste(myOutDir_3three_diffMe,  "1B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_3three_sub3_1kb, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_3three_sub3_1kb, nclass=100, xlim=c(0, 50), freq=FALSE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
hist(myMethDi_3three_sub3_1kb, nclass=100, xlim=c(0, 100), freq=TRUE)
hist(myMethDi_3three_sub3_1kb, nclass=100, xlim=c(0, 50), freq=TRUE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb>=20], nclass=81, xlim=c(20, 100),  freq=TRUE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb>=10], nclass=91, xlim=c(10, 100),  freq=TRUE)
hist(myMethDi_3three_sub3_1kb[myMethDi_3three_sub3_1kb<=10], nclass=100, xlim=c(0, 10),   freq=TRUE)
dev.off() 




myColor1_3three_sub3_1kb <- rep( "no",   times= length(myQvalue_3three_sub3_1kb) )
myColor1_3three_sub3_1kb[ (abs(myMethDi_3three_sub3_1kb)>10) & (myQvalue_3three_sub3_1kb<0.001) ]  <- "yes"
length( myColor1_3three_sub3_1kb[myColor1_3three_sub3_1kb=="yes"] )



DataFrame2_3three_sub3_1kb <- data.frame(myx1 =  myMethDi_3three_sub3_1kb,   
                                         myy1 =  -log10(myQvalue_3three_sub3_1kb),  
                                         mycolor1 = myColor1_3three_sub3_1kb )

ggplot(DataFrame2_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_diffMe,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_3three_diffMe,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_3three_diffMe,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )




ggplot(DataFrame2_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_3three_diffMe,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_3three_diffMe,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_3three_sub3_1kb, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_3three_diffMe,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )









myDiff25p.hypo_3three_sub3_1kb  = getMethylDiff(myDiff_3three_sub3_1kb, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_3three_sub3_1kb = getMethylDiff(myDiff_3three_sub3_1kb, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
myDiff25p_3three_sub3_1kb       = getMethylDiff(myDiff_3three_sub3_1kb, difference=10, qvalue=0.05)
myDiffTemp_3three_sub3_1kb      = getMethylDiff(myDiff_3three_sub3_1kb, difference=0,  qvalue=0.05)

write.table(myDiff_3three_sub3_1kb , file = paste(myOutDir_3three_diffMe,"2A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_3three_sub3_1kb , file = paste(myOutDir_3three_diffMe,"2B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_3three_sub3_1kb , file = paste(myOutDir_3three_diffMe,"2C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_3three_sub3_1kb , file = paste(myOutDir_3three_diffMe,"2D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_3three_sub3_1kb , file = paste(myOutDir_3three_diffMe,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_3three_diffMe, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_3three_sub3_1kb,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_3three_diffMe, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_3three_sub3_1kb,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
diffGeneAnn_3three_sub3_1kb_hypo = annotateWithGeneParts(as(myDiff25p.hypo_3three_sub3_1kb,  "GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_diffMe, "4A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffGeneAnn_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub3_1kb_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub3_1kb,"GRanges"),
                                                           cpg.obj$CpGi,  cpg.obj$shores,
                                                           feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "4C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffCpGann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_3three_sub3_1kb,"GRanges"),
                                                              myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                              feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "4E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffrepeatann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "5A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "5C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "5E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "5G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub3_1kb_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_3three_sub3_1kb, "GRanges"),
                                                                  feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                  feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "5I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub3_1kb_hypo,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub3_1kb_hypo)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub3_1kb_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_3three_sub3_1kb_hyper = annotateWithGeneParts(as(myDiff25p.hyper_3three_sub3_1kb,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_3three_diffMe, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffGeneAnn_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_3three_sub3_1kb_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub3_1kb,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffCpGann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_3three_sub3_1kb,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffrepeatann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted1Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted2Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted3Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted4Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_3three_sub3_1kb_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_3three_sub3_1kb, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_3three_diffMe, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_3three_sub3_1kb_hyper,  percentage=TRUE)
print(diffImprinted5Ann_3three_sub3_1kb_hyper)
sink()

pdf( file=paste(myOutDir_3three_diffMe, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_3three_sub3_1kb_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##################################################################################################################
##################################################################################################################








































######################################################################################################################################################
######################################################################################################################################################
myOutDir_4four = paste(myOutDir, "/4-500bpBin",  sep="");
if( ! file.exists(myOutDir_4four) ) { dir.create(myOutDir_4four, recursive = TRUE) }


#####
tiles_4four_sub4_500bp = tileMethylCounts( myobj, win.size=500, step.size=500) ## 500bp bin
meth_4four_sub4_500bp  = unite( tiles_4four_sub4_500bp, destrand=FALSE   )   ## 100% overlap
mat_4four_sub4_500bp   = percMethylation( meth_4four_sub4_500bp )


sink( file=paste(myOutDir_4four , "0-dimensions-500bpBin.txt", sep="/")  )
print( tiles_4four_sub4_500bp )
print("#########dimensions:")
print( dim(meth_4four_sub4_500bp)  )   
print( dim(mat_4four_sub4_500bp)   )
sink()




write.table(meth_4four_sub4_500bp , 
            file = paste(myOutDir_4four,   "1A-meth-500bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_4four_sub4_500bp , 
            file = paste(myOutDir_4four,   "1B-mat-500bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





pdf( file=paste(myOutDir_4four, "2A-MethylationStats-500bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_4four_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_4four, "2B-MethylationStats-500bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_4four_sub4_500bp[[i]] )  )
}
sink()




pdf( file=paste(myOutDir_4four, "3A-CoverageStats-500bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_4four_sub4_500bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_4four, "3B-CoverageStats-500bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_4four_sub4_500bp[[i]] )  )
}
sink()





MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4A-cluster-kept100percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0 )

MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4B-cluster-kept80percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.2 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4C-cluster-kept60percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.4 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4D-cluster-kept50percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.5 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4E-cluster-kept40percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.6 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4F-cluster-kept30percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.7 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4G-cluster-kept20percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.8 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4H-cluster-kept10percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.9 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4I-cluster-kept5percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.95 )


MyCluster_1(  mymeth1 = meth_4four_sub4_500bp , path1 = myOutDir_4four, 
              file1 = "4J-cluster-kept1percent-500bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.99 )





pdf( file=paste(myOutDir_4four, "5-PCA-500bpBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_4four_sub4_500bp , screeplot=TRUE)
PCASamples(meth_4four_sub4_500bp )
dev.off()

#####################





PCA_4four_sub4_500bp <- prcomp( t(mat_4four_sub4_500bp)  )
names(PCA_4four_sub4_500bp)

sink( file = paste(myOutDir_4four,  "6A-PCA-500bpBin.txt",  sep="/") )
print(PCA_4four_sub4_500bp)
sink()

sink( file = paste(myOutDir_4four,  "6B-PCA-summary-500bpBin.txt",  sep="/") )
summary(PCA_4four_sub4_500bp)
sink()


sink( file = paste(myOutDir_4four,  "6C-PCA-all-500bpBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_4four_sub4_500bp$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_4four_sub4_500bp$rotation)
print("####################### myPCA2$center #########################")
print(PCA_4four_sub4_500bp$center)
print("####################### myPCA2$scale #########################")
print(PCA_4four_sub4_500bp$scale)
print("####################### myPCA2$x #########################")
print(PCA_4four_sub4_500bp$x)
sink()


pdf( file=paste(myOutDir_4four,   "6D-PCA-info-500bpBin.pdf", sep="/")  )
plot(PCA_4four_sub4_500bp, type="lines")
fviz_eig(PCA_4four_sub4_500bp)
dev.off() 




my_fviz_pca_ind1_4four_sub4_500bp <- fviz_pca_ind(PCA_4four_sub4_500bp,
                                                  col.ind = "cos2", # Color by the quality of representation
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_4four_sub4_500bp <- fviz_pca_ind(PCA_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE
)
my_fviz_pca_ind3_4four_sub4_500bp <- fviz_pca_ind(PCA_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ) , # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)
my_fviz_pca_ind4_4four_sub4_500bp <- fviz_pca_ind(PCA_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  #legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)


svg(file=paste(myOutDir_4four, "7A-PCA-2D-1-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "7B-PCA-2D-2-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "7C-PCA-2D-3-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "7D-PCA-2D-4-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_4four_sub4_500bp)
dev.off() 








#############################


PCA_4four_sub4_500bp_matrix <- PCA_4four_sub4_500bp$x
dim( PCA_4four_sub4_500bp_matrix )

PCA_4four_sub4_500bp_Contri  <- (PCA_4four_sub4_500bp$sdev)^2
PCA_4four_sub4_500bp_Contri  <- PCA_4four_sub4_500bp_Contri/sum(PCA_4four_sub4_500bp_Contri)
PCA_4four_sub4_500bp_Contri  <- PCA_4four_sub4_500bp_Contri * 100
PCA_4four_sub4_500bp_Contri  <- round(PCA_4four_sub4_500bp_Contri, 2)

label1_4four_sub4_500bp <-   paste( "PC1 ",  "(", PCA_4four_sub4_500bp_Contri[1], "%)", sep="" )
label2_4four_sub4_500bp <-   paste( "PC2 ",  "(", PCA_4four_sub4_500bp_Contri[2], "%)", sep="" )
label3_4four_sub4_500bp <-   paste( "PC3 ",  "(", PCA_4four_sub4_500bp_Contri[3], "%)", sep="" )
label1_4four_sub4_500bp  
label2_4four_sub4_500bp  
label3_4four_sub4_500bp  


myLabel = myType1
dataframeA_4four_sub4_500bp  <- data.frame( as.data.frame(PCA_4four_sub4_500bp_matrix), myType2, myType3, myLabel   ) 
dataframeA_4four_sub4_500bp 


FigureTemp1_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="8A-PCA-PC1-PC2-500bpBin",  height1=3,  width1=4.8)


FigureTemp2_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="8B-PCA-PC1-PC2-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp3_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="8C-PCA-PC1-PC2-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="8D-PCA-PC1-PC2-big-500bpBin",   height1=3,  width1=4.8)



FigureTemp5_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="8E-PCA-PC1-PC2-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp6_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="8F-PCA-PC1-PC2-text2-500bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="9A-PCA-PC1-PC3-500bpBin",   height1=3,  width1=4.8)


FigureTemp12_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="9B-PCA-PC1-PC3-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp13_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="9C-PCA-PC1-PC3-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="9D-PCA-PC1-PC3-big-500bpBin",  height1=3,  width1=4.8)



FigureTemp15_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="9E-PCA-PC1-PC3-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp16_4four_sub4_500bp  <- ggplot( data = dataframeA_4four_sub4_500bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="9F-PCA-PC1-PC3-text2-500bpBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_4four, "10A_PCA-3d-500bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_4four_sub4_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_4four_sub4_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_4four_sub4_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_4four, "10B_PCA-3d-500bpBin-withLabels.pdf",  sep="/") )

s3d1_4four_sub4_500bp <- scatterplot3d( 
  dataframeA_4four_sub4_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_4four_sub4_500bp$xyz.convert(dataframeA_4four_sub4_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_4four_sub4_500bp <- scatterplot3d( 
  dataframeA_4four_sub4_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_4four_sub4_500bp$xyz.convert(dataframeA_4four_sub4_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_4four_sub4_500bp <- scatterplot3d( 
  dataframeA_4four_sub4_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_4four_sub4_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_4four_sub4_500bp$xyz.convert(dataframeA_4four_sub4_500bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "11A_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "11B_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_4four, "12_PCA-3d-500bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )



##
scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )



##
scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeA_4four_sub4_500bp[,1], y = dataframeA_4four_sub4_500bp[,2], z = dataframeA_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.


myres.pca_4four_sub4_500bp  <- PCA(t(mat_4four_sub4_500bp),   graph = FALSE)

sink( file = paste(myOutDir_4four, "20_PCA-info-500bpBin-byPCA.txt",  sep="/") )
print( myres.pca_4four_sub4_500bp )
print( "#################################" )
print( summary(myres.pca_4four_sub4_500bp) )
print( "#################################"  )
myeig.val_4four_sub4_500bp <- get_eigenvalue( myres.pca_4four_sub4_500bp )
myeig.val_4four_sub4_500bp
sink() 


pdf( file = paste(myOutDir_4four, "21_PCA-screePlot-500bpBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_4four_sub4_500bp, addlabels = TRUE )
fviz_screeplot(X=myres.pca_4four_sub4_500bp, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_4four_sub4_500bp, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_4four_sub4_500bp, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_4four_sub4_500bp, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_4four_sub4_500bp <- fviz_pca_ind(myres.pca_4four_sub4_500bp,
                                                  col.ind = "cos2", # Color by the quality of representation
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_4four_sub4_500bp <- fviz_pca_ind(myres.pca_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE
)
my_fviz_pca_ind3_4four_sub4_500bp <- fviz_pca_ind(myres.pca_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ) , # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)
my_fviz_pca_ind4_4four_sub4_500bp <- fviz_pca_ind(myres.pca_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  #legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)
my_fviz_pca_ind5_4four_sub4_500bp <- fviz_pca_ind(myres.pca_4four_sub4_500bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE
)

svg(file=paste(myOutDir_4four, "22A-PCA-2D-1-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "22B-PCA-2D-2-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "22C-PCA-2D-3-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "22D-PCA-2D-4-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_4four_sub4_500bp)
dev.off() 

svg(file=paste(myOutDir_4four, "22E-PCA-2D-4-500bpBin.svg", sep="/") )
print(my_fviz_pca_ind5_4four_sub4_500bp)
dev.off() 



#############################


myres.pca_4four_sub4_500bp_matrix <- myres.pca_4four_sub4_500bp$ind$coord
dim( myres.pca_4four_sub4_500bp_matrix )

myres.pca_4four_sub4_500bp_Contri  <- (myres.pca_4four_sub4_500bp$eig)[,2] 
myres.pca_4four_sub4_500bp_Contri  <- myres.pca_4four_sub4_500bp_Contri/sum(myres.pca_4four_sub4_500bp_Contri)
myres.pca_4four_sub4_500bp_Contri  <- myres.pca_4four_sub4_500bp_Contri * 100
myres.pca_4four_sub4_500bp_Contri  <- round(myres.pca_4four_sub4_500bp_Contri, 2)

label1_4four_sub4_500bp <-   paste( "PC1 ",  "(", myres.pca_4four_sub4_500bp_Contri[1], "%)", sep="" )
label2_4four_sub4_500bp <-   paste( "PC2 ",  "(", myres.pca_4four_sub4_500bp_Contri[2], "%)", sep="" )
label3_4four_sub4_500bp <-   paste( "PC3 ",  "(", myres.pca_4four_sub4_500bp_Contri[3], "%)", sep="" )
label1_4four_sub4_500bp  
label2_4four_sub4_500bp  
label3_4four_sub4_500bp  


myLabel = myType1
dataframeB_4four_sub4_500bp  <- data.frame( as.data.frame(myres.pca_4four_sub4_500bp_matrix), myType2, myType3, myLabel   ) 
dataframeB_4four_sub4_500bp 


FigureTemp1_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="23A-PCA-PC1-PC2-500bpBin",  height1=3,  width1=4.8)


FigureTemp2_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="23B-PCA-PC1-PC2-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp3_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="23C-PCA-PC1-PC2-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="23D-PCA-PC1-PC2-big-500bpBin",   height1=3,  width1=4.8)



FigureTemp5_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="23E-PCA-PC1-PC2-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp6_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label2_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="23F-PCA-PC1-PC2-text2-500bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="24A-PCA-PC1-PC3-500bpBin",   height1=3,  width1=4.8)


FigureTemp12_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="24B-PCA-PC1-PC3-alpha-500bpBin",   height1=3,  width1=4.8)


FigureTemp13_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="24C-PCA-PC1-PC3-smallDot-500bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="24D-PCA-PC1-PC3-big-500bpBin",  height1=3,  width1=4.8)



FigureTemp15_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="24E-PCA-PC1-PC3-text-500bpBin",   height1=3,  width1=4.8)


FigureTemp16_4four_sub4_500bp  <- ggplot( data = dataframeB_4four_sub4_500bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_4four_sub4_500bp) +   ylab(label3_4four_sub4_500bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_4four_sub4_500bp ,  path1=myOutDir_4four, fileName1="24F-PCA-PC1-PC3-text2-500bpBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_4four, "25A_PCA-3d-500bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_4four_sub4_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_4four_sub4_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_4four_sub4_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_4four, "25B_PCA-3d-500bpBin-withLabels.pdf",  sep="/") )

s3d1_4four_sub4_500bp <- scatterplot3d( 
  dataframeB_4four_sub4_500bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_500bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_4four_sub4_500bp$xyz.convert(dataframeB_4four_sub4_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_4four_sub4_500bp <- scatterplot3d( 
  dataframeB_4four_sub4_500bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_500bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_4four_sub4_500bp$xyz.convert(dataframeB_4four_sub4_500bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_4four_sub4_500bp <- scatterplot3d( 
  dataframeB_4four_sub4_500bp[,1:3], pch =myType2_shape, 
  xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_4four_sub4_500bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_4four_sub4_500bp$xyz.convert(dataframeB_4four_sub4_500bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "26A_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_4four, "26B_PCA-3d-500bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_4four, "27_PCA-3d-500bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )



##
scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )



##
scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )


scatter3D(x = dataframeB_4four_sub4_500bp[,1], y = dataframeB_4four_sub4_500bp[,2], z = dataframeB_4four_sub4_500bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_4four_sub4_500bp,  ylab = label2_4four_sub4_500bp,  zlab = label3_4four_sub4_500bp  )

dev.off()






#################



km.res1_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=2, iter.max = 20 )
km.res2_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=3, iter.max = 20  )
km.res3_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=4, iter.max = 20  )
km.res4_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=5, iter.max = 20 )
km.res5_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=6, iter.max = 20 )
km.res6_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=7, iter.max = 20 )
km.res7_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=8, iter.max = 20 )
km.res8_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=9, iter.max = 20 )
km.res9_4four_sub4_500bp  <- kmeans( t(mat_4four_sub4_500bp), centers=10, iter.max = 20 )
km.res10_4four_sub4_500bp <- kmeans( t(mat_4four_sub4_500bp), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_4four, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_4four_sub4_500bp$cluster
print("#################### 3 classes:")
km.res2_4four_sub4_500bp$cluster
print("#################### 4 classes:")
km.res3_4four_sub4_500bp$cluster
print("#################### 5 classes:")
km.res4_4four_sub4_500bp$cluster
print("#################### 6 classes:")
km.res5_4four_sub4_500bp$cluster
print("#################### 7 classes:")
km.res6_4four_sub4_500bp$cluster
print("#################### 8 classes:")
km.res7_4four_sub4_500bp$cluster
print("#################### 9 classes:")
km.res8_4four_sub4_500bp$cluster
print("#################### 10 classes:")
km.res9_4four_sub4_500bp$cluster
print("#################### 11 classes:")
km.res10_4four_sub4_500bp$cluster
sink()



pdf( file = paste(myOutDir_4four, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_500bp), col = km.res1_4four_sub4_500bp$cluster )

plot( t(mat_4four_sub4_500bp), col = km.res1_4four_sub4_500bp$cluster)
points(km.res1_4four_sub4_500bp$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_500bp), col = km.res1_4four_sub4_500bp$cluster)
points(km.res1_4four_sub4_500bp$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_4four_sub4_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 






pdf( file = paste(myOutDir_4four, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_500bp), col = km.res2_4four_sub4_500bp$cluster )

plot( t(mat_4four_sub4_500bp), col = km.res2_4four_sub4_500bp$cluster)
points(km.res2_4four_sub4_500bp$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_500bp), col = km.res2_4four_sub4_500bp$cluster)
points(km.res2_4four_sub4_500bp$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_4four_sub4_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_4four, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_500bp), col = km.res3_4four_sub4_500bp$cluster )

plot( t(mat_4four_sub4_500bp), col = km.res3_4four_sub4_500bp$cluster)
points(km.res3_4four_sub4_500bp$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_500bp), col = km.res3_4four_sub4_500bp$cluster)
points(km.res3_4four_sub4_500bp$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_4four_sub4_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 







pdf( file = paste(myOutDir_4four, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_4four_sub4_500bp), col = km.res4_4four_sub4_500bp$cluster )

plot( t(mat_4four_sub4_500bp), col = km.res4_4four_sub4_500bp$cluster)
points(km.res4_4four_sub4_500bp$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_4four_sub4_500bp), col = km.res4_4four_sub4_500bp$cluster)
points(km.res4_4four_sub4_500bp$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_4four_sub4_500bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "euclidean")
res.dist2_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "maximum")
res.dist3_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "manhattan")
res.dist4_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "canberra")
res.dist5_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "binary")
res.dist6_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "minkowski")
res.dist7_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "pearson")
res.dist8_4four_sub4_500bp <- get_dist( t(mat_4four_sub4_500bp) ,   method = "spearman")




pdf( file = paste(myOutDir_4four, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_4four_sub4_500bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_4four, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_4four_sub4_500bp ) 
print("################### maximum: ")
print(res.dist2_4four_sub4_500bp ) 
print("################### manhattan: ")
print(res.dist3_4four_sub4_500bp ) 
print("################### canberra: ")
print(res.dist4_4four_sub4_500bp ) 
print("################### binary: ")
print(res.dist5_4four_sub4_500bp ) 
print("################### minkowski: ")
print(res.dist6_4four_sub4_500bp ) 
print("################### pearson: ")
print(res.dist7_4four_sub4_500bp ) 
print("################### spearman: ")
print(res.dist8_4four_sub4_500bp ) 
sink() 


# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "ward.D" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "ward.D" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "ward.D" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "ward.D" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "ward.D" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "ward.D" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "ward.D" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "ward.D" )  

pdf( file = paste(myOutDir_4four, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "ward.D2" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "ward.D2" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "ward.D2" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "ward.D2" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "ward.D2" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "ward.D2" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "ward.D2" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "ward.D2" )  

pdf( file = paste(myOutDir_4four, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "single" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "single" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "single" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "single" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "single" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "single" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "single" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "single" )  

pdf( file = paste(myOutDir_4four, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "complete" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "complete" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "complete" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "complete" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "complete" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "complete" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "complete" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "complete" )  

pdf( file = paste(myOutDir_4four, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "average" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "average" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "average" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "average" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "average" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "average" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "average" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "average" )  

pdf( file = paste(myOutDir_4four, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "mcquitty" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "mcquitty" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "mcquitty" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "mcquitty" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "mcquitty" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "mcquitty" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "mcquitty" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "mcquitty" )  

pdf( file = paste(myOutDir_4four, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "median" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "median" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "median" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "median" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "median" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "median" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "median" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "median" )  

pdf( file = paste(myOutDir_4four, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_4four_sub4_500bp <- hclust(res.dist1_4four_sub4_500bp, method = "centroid" )   
res.hc2_4four_sub4_500bp <- hclust(res.dist2_4four_sub4_500bp, method = "centroid" )   
res.hc3_4four_sub4_500bp <- hclust(res.dist3_4four_sub4_500bp, method = "centroid" )   
res.hc4_4four_sub4_500bp <- hclust(res.dist4_4four_sub4_500bp, method = "centroid" )   
res.hc5_4four_sub4_500bp <- hclust(res.dist5_4four_sub4_500bp, method = "centroid" )   
res.hc6_4four_sub4_500bp <- hclust(res.dist6_4four_sub4_500bp, method = "centroid" )   
res.hc7_4four_sub4_500bp <- hclust(res.dist7_4four_sub4_500bp, method = "centroid" )   
res.hc8_4four_sub4_500bp <- hclust(res.dist8_4four_sub4_500bp, method = "centroid" )  

pdf( file = paste(myOutDir_4four, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_4four_sub4_500bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_4four, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp )
fviz_dend(res.hc2_4four_sub4_500bp )
fviz_dend(res.hc3_4four_sub4_500bp )
fviz_dend(res.hc4_4four_sub4_500bp )
fviz_dend(res.hc5_4four_sub4_500bp )
fviz_dend(res.hc6_4four_sub4_500bp )
fviz_dend(res.hc7_4four_sub4_500bp )
fviz_dend(res.hc8_4four_sub4_500bp )
dev.off() 


pdf( file = paste(myOutDir_4four, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc2_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc3_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc4_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc5_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc6_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc7_4four_sub4_500bp, label_cols = myType3_color )
fviz_dend(res.hc8_4four_sub4_500bp, label_cols = myType3_color )
dev.off() 









##################
myOutDir_4four_diffMe   <- paste(myOutDir_4four, "diffMe_DMC_DMR",  sep="/")
if( ! file.exists(myOutDir_4four_diffMe) ) { dir.create(myOutDir_4four_diffMe, recursive = TRUE) }


myDiff_4four_sub4_500bp = calculateDiffMeth(meth_4four_sub4_500bp,  num.cores=16)
dim(myDiff_4four_sub4_500bp)
names(myDiff_4four_sub4_500bp)
head(myDiff_4four_sub4_500bp)

myQvalue_4four_sub4_500bp = myDiff_4four_sub4_500bp$qvalue
myMethDi_4four_sub4_500bp = myDiff_4four_sub4_500bp$meth.diff



pdf(paste(myOutDir_4four_diffMe,  "1A_qvalue_distribution.pdf",  sep="/"))
hist( myQvalue_4four_sub4_500bp,  nclass=100, xlim=c(0, 1), freq=FALSE)
hist( myQvalue_4four_sub4_500bp[myQvalue_4four_sub4_500bp<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE)
hist( myQvalue_4four_sub4_500bp,  nclass=100, xlim=c(0, 1), freq=TRUE)
hist( myQvalue_4four_sub4_500bp[myQvalue_4four_sub4_500bp<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE)
dev.off() 


pdf(paste(myOutDir_4four_diffMe,  "1B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_4four_sub4_500bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_4four_sub4_500bp, nclass=100, xlim=c(0, 50), freq=FALSE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
hist(myMethDi_4four_sub4_500bp, nclass=100, xlim=c(0, 100), freq=TRUE)
hist(myMethDi_4four_sub4_500bp, nclass=100, xlim=c(0, 50), freq=TRUE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp>=20], nclass=81, xlim=c(20, 100),  freq=TRUE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp>=10], nclass=91, xlim=c(10, 100),  freq=TRUE)
hist(myMethDi_4four_sub4_500bp[myMethDi_4four_sub4_500bp<=10], nclass=100, xlim=c(0, 10),   freq=TRUE)
dev.off() 




myColor1_4four_sub4_500bp <- rep( "no",   times= length(myQvalue_4four_sub4_500bp) )
myColor1_4four_sub4_500bp[ (abs(myMethDi_4four_sub4_500bp)>10) & (myQvalue_4four_sub4_500bp<0.001) ]  <- "yes"
length( myColor1_4four_sub4_500bp[myColor1_4four_sub4_500bp=="yes"] )



DataFrame2_4four_sub4_500bp <- data.frame(myx1 =  myMethDi_4four_sub4_500bp,   
                                          myy1 =  -log10(myQvalue_4four_sub4_500bp),  
                                          mycolor1 = myColor1_4four_sub4_500bp )

ggplot(DataFrame2_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_diffMe,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_4four_diffMe,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_4four_diffMe,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )




ggplot(DataFrame2_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_4four_diffMe,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_4four_diffMe,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_4four_sub4_500bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_4four_diffMe,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )









myDiff25p.hypo_4four_sub4_500bp  = getMethylDiff(myDiff_4four_sub4_500bp, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_4four_sub4_500bp = getMethylDiff(myDiff_4four_sub4_500bp, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
myDiff25p_4four_sub4_500bp       = getMethylDiff(myDiff_4four_sub4_500bp, difference=10, qvalue=0.05)
myDiffTemp_4four_sub4_500bp      = getMethylDiff(myDiff_4four_sub4_500bp, difference=0,  qvalue=0.05)

write.table(myDiff_4four_sub4_500bp , file = paste(myOutDir_4four_diffMe,"2A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_4four_sub4_500bp , file = paste(myOutDir_4four_diffMe,"2B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_4four_sub4_500bp , file = paste(myOutDir_4four_diffMe,"2C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_4four_sub4_500bp , file = paste(myOutDir_4four_diffMe,"2D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_4four_sub4_500bp , file = paste(myOutDir_4four_diffMe,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_4four_diffMe, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_4four_sub4_500bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_4four_diffMe, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_4four_sub4_500bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
diffGeneAnn_4four_sub4_500bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_4four_sub4_500bp,  "GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_diffMe, "4A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffGeneAnn_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub4_500bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub4_500bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "4C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffCpGann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_4four_sub4_500bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "4E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffrepeatann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "5A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "5C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "5E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "5G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub4_500bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_4four_sub4_500bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "5I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub4_500bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub4_500bp_hypo)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub4_500bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_4four_sub4_500bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_4four_sub4_500bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_4four_diffMe, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffGeneAnn_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_4four_sub4_500bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub4_500bp,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffCpGann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_4four_sub4_500bp,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffrepeatann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_4four_sub4_500bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_4four_sub4_500bp, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_4four_diffMe, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_4four_sub4_500bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_4four_sub4_500bp_hyper)
sink()

pdf( file=paste(myOutDir_4four_diffMe, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_4four_sub4_500bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##################################################################################################################
##################################################################################################################





































######################################################################################################################################################
######################################################################################################################################################
myOutDir_5five = paste(myOutDir, "/5-100bpBin",  sep="");
if( ! file.exists(myOutDir_5five) ) { dir.create(myOutDir_5five, recursive = TRUE) }


#####
tiles_5five_sub5_100bp = tileMethylCounts( myobj, win.size=100, step.size=100) ## 100bp bin
meth_5five_sub5_100bp  = unite( tiles_5five_sub5_100bp, destrand=FALSE   )   ## 100% overlap
mat_5five_sub5_100bp   = percMethylation( meth_5five_sub5_100bp )


sink( file=paste(myOutDir_5five , "0-dimensions-100bpBin.txt", sep="/")  )
print( tiles_5five_sub5_100bp )
print("#########dimensions:")
print( dim(meth_5five_sub5_100bp)  )   
print( dim(mat_5five_sub5_100bp)   )
sink()




write.table(meth_5five_sub5_100bp , 
            file = paste(myOutDir_5five,   "1A-meth-100bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_5five_sub5_100bp , 
            file = paste(myOutDir_5five,   "1B-mat-100bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





pdf( file=paste(myOutDir_5five, "2A-MethylationStats-100bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_5five_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_5five, "2B-MethylationStats-100bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_5five_sub5_100bp[[i]] )  )
}
sink()




pdf( file=paste(myOutDir_5five, "3A-CoverageStats-100bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_5five_sub5_100bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_5five, "3B-CoverageStats-100bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_5five_sub5_100bp[[i]] )  )
}
sink()





MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4A-cluster-kept100percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0 )

MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4B-cluster-kept80percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.2 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4C-cluster-kept60percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.4 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4D-cluster-kept50percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.5 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4E-cluster-kept40percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.6 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4F-cluster-kept30percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.7 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4G-cluster-kept20percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.8 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4H-cluster-kept10percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.9 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4I-cluster-kept5percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.95 )


MyCluster_1(  mymeth1 = meth_5five_sub5_100bp , path1 = myOutDir_5five, 
              file1 = "4J-cluster-kept1percent-100bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.99 )





pdf( file=paste(myOutDir_5five, "5-PCA-100bpBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_5five_sub5_100bp , screeplot=TRUE)
PCASamples(meth_5five_sub5_100bp )
dev.off()

#####################





PCA_5five_sub5_100bp <- prcomp( t(mat_5five_sub5_100bp)  )
names(PCA_5five_sub5_100bp)

sink( file = paste(myOutDir_5five,  "6A-PCA-100bpBin.txt",  sep="/") )
print(PCA_5five_sub5_100bp)
sink()

sink( file = paste(myOutDir_5five,  "6B-PCA-summary-100bpBin.txt",  sep="/") )
summary(PCA_5five_sub5_100bp)
sink()


sink( file = paste(myOutDir_5five,  "6C-PCA-all-100bpBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_5five_sub5_100bp$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_5five_sub5_100bp$rotation)
print("####################### myPCA2$center #########################")
print(PCA_5five_sub5_100bp$center)
print("####################### myPCA2$scale #########################")
print(PCA_5five_sub5_100bp$scale)
print("####################### myPCA2$x #########################")
print(PCA_5five_sub5_100bp$x)
sink()


pdf( file=paste(myOutDir_5five,   "6D-PCA-info-100bpBin.pdf", sep="/")  )
plot(PCA_5five_sub5_100bp, type="lines")
fviz_eig(PCA_5five_sub5_100bp)
dev.off() 




my_fviz_pca_ind1_5five_sub5_100bp <- fviz_pca_ind(PCA_5five_sub5_100bp,
                                                  col.ind = "cos2", # Color by the quality of representation
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_5five_sub5_100bp <- fviz_pca_ind(PCA_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE
)
my_fviz_pca_ind3_5five_sub5_100bp <- fviz_pca_ind(PCA_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ) , # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)
my_fviz_pca_ind4_5five_sub5_100bp <- fviz_pca_ind(PCA_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  #legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)


svg(file=paste(myOutDir_5five, "7A-PCA-2D-1-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "7B-PCA-2D-2-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "7C-PCA-2D-3-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "7D-PCA-2D-4-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_5five_sub5_100bp)
dev.off() 








#############################


PCA_5five_sub5_100bp_matrix <- PCA_5five_sub5_100bp$x
dim( PCA_5five_sub5_100bp_matrix )

PCA_5five_sub5_100bp_Contri  <- (PCA_5five_sub5_100bp$sdev)^2
PCA_5five_sub5_100bp_Contri  <- PCA_5five_sub5_100bp_Contri/sum(PCA_5five_sub5_100bp_Contri)
PCA_5five_sub5_100bp_Contri  <- PCA_5five_sub5_100bp_Contri * 100
PCA_5five_sub5_100bp_Contri  <- round(PCA_5five_sub5_100bp_Contri, 2)

label1_5five_sub5_100bp <-   paste( "PC1 ",  "(", PCA_5five_sub5_100bp_Contri[1], "%)", sep="" )
label2_5five_sub5_100bp <-   paste( "PC2 ",  "(", PCA_5five_sub5_100bp_Contri[2], "%)", sep="" )
label3_5five_sub5_100bp <-   paste( "PC3 ",  "(", PCA_5five_sub5_100bp_Contri[3], "%)", sep="" )
label1_5five_sub5_100bp  
label2_5five_sub5_100bp  
label3_5five_sub5_100bp  


myLabel = myType1
dataframeA_5five_sub5_100bp  <- data.frame( as.data.frame(PCA_5five_sub5_100bp_matrix), myType2, myType3, myLabel   ) 
dataframeA_5five_sub5_100bp 


FigureTemp1_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="8A-PCA-PC1-PC2-100bpBin",  height1=3,  width1=4.8)


FigureTemp2_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="8B-PCA-PC1-PC2-alpha-100bpBin",   height1=3,  width1=4.8)


FigureTemp3_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="8C-PCA-PC1-PC2-smallDot-100bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="8D-PCA-PC1-PC2-big-100bpBin",   height1=3,  width1=4.8)



FigureTemp5_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="8E-PCA-PC1-PC2-text-100bpBin",   height1=3,  width1=4.8)


FigureTemp6_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="8F-PCA-PC1-PC2-text2-100bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="9A-PCA-PC1-PC3-100bpBin",   height1=3,  width1=4.8)


FigureTemp12_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="9B-PCA-PC1-PC3-alpha-100bpBin",   height1=3,  width1=4.8)


FigureTemp13_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="9C-PCA-PC1-PC3-smallDot-100bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="9D-PCA-PC1-PC3-big-100bpBin",  height1=3,  width1=4.8)



FigureTemp15_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="9E-PCA-PC1-PC3-text-100bpBin",   height1=3,  width1=4.8)


FigureTemp16_5five_sub5_100bp  <- ggplot( data = dataframeA_5five_sub5_100bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="9F-PCA-PC1-PC3-text2-100bpBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_5five, "10A_PCA-3d-100bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_5five_sub5_100bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_5five_sub5_100bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_5five_sub5_100bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_5five_sub5_100bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_5five_sub5_100bp[,1:3], pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_5five_sub5_100bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_5five, "10B_PCA-3d-100bpBin-withLabels.pdf",  sep="/") )

s3d1_5five_sub5_100bp <- scatterplot3d( 
  dataframeA_5five_sub5_100bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_5five_sub5_100bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_5five_sub5_100bp$xyz.convert(dataframeA_5five_sub5_100bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_5five_sub5_100bp <- scatterplot3d( 
  dataframeA_5five_sub5_100bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_5five_sub5_100bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_5five_sub5_100bp$xyz.convert(dataframeA_5five_sub5_100bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_5five_sub5_100bp <- scatterplot3d( 
  dataframeA_5five_sub5_100bp[,1:3], pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_5five_sub5_100bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_5five_sub5_100bp$xyz.convert(dataframeA_5five_sub5_100bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_5five, "11A_PCA-3d-100bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_5five, "11B_PCA-3d-100bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_5five, "12_PCA-3d-100bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )



##
scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )



##
scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeA_5five_sub5_100bp[,1], y = dataframeA_5five_sub5_100bp[,2], z = dataframeA_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.


myres.pca_5five_sub5_100bp  <- PCA(t(mat_5five_sub5_100bp),   graph = FALSE)

sink( file = paste(myOutDir_5five, "20_PCA-info-100bpBin-byPCA.txt",  sep="/") )
print( myres.pca_5five_sub5_100bp )
print( "#################################" )
print( summary(myres.pca_5five_sub5_100bp) )
print( "#################################"  )
myeig.val_5five_sub5_100bp <- get_eigenvalue( myres.pca_5five_sub5_100bp )
myeig.val_5five_sub5_100bp
sink() 


pdf( file = paste(myOutDir_5five, "21_PCA-screePlot-100bpBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_5five_sub5_100bp, addlabels = TRUE )
fviz_screeplot(X=myres.pca_5five_sub5_100bp, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_5five_sub5_100bp, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_5five_sub5_100bp, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_5five_sub5_100bp, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_5five_sub5_100bp <- fviz_pca_ind(myres.pca_5five_sub5_100bp,
                                                  col.ind = "cos2", # Color by the quality of representation
                                                  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                                  repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_5five_sub5_100bp <- fviz_pca_ind(myres.pca_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE
)
my_fviz_pca_ind3_5five_sub5_100bp <- fviz_pca_ind(myres.pca_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ) , # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  addEllipses = TRUE, # Concentration ellipses
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)
my_fviz_pca_ind4_5five_sub5_100bp <- fviz_pca_ind(myres.pca_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  #legend.title = "Groups",
                                                  repel = TRUE, 
                                                  label = "none", 
                                                  alpha.ind = 1
)
my_fviz_pca_ind5_5five_sub5_100bp <- fviz_pca_ind(myres.pca_5five_sub5_100bp,
                                                  col.ind =  as.factor( myTreatment ), # color by groups
                                                  #palette = c("#00AFBB",  "#FC4E07"),
                                                  ellipse.type = "confidence",
                                                  legend.title = "Groups",
                                                  repel = TRUE
)

svg(file=paste(myOutDir_5five, "22A-PCA-2D-1-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "22B-PCA-2D-2-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "22C-PCA-2D-3-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "22D-PCA-2D-4-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_5five_sub5_100bp)
dev.off() 

svg(file=paste(myOutDir_5five, "22E-PCA-2D-4-100bpBin.svg", sep="/") )
print(my_fviz_pca_ind5_5five_sub5_100bp)
dev.off() 



#############################


myres.pca_5five_sub5_100bp_matrix <- myres.pca_5five_sub5_100bp$ind$coord
dim( myres.pca_5five_sub5_100bp_matrix )

myres.pca_5five_sub5_100bp_Contri  <- (myres.pca_5five_sub5_100bp$eig)[,2] 
myres.pca_5five_sub5_100bp_Contri  <- myres.pca_5five_sub5_100bp_Contri/sum(myres.pca_5five_sub5_100bp_Contri)
myres.pca_5five_sub5_100bp_Contri  <- myres.pca_5five_sub5_100bp_Contri * 100
myres.pca_5five_sub5_100bp_Contri  <- round(myres.pca_5five_sub5_100bp_Contri, 2)

label1_5five_sub5_100bp <-   paste( "PC1 ",  "(", myres.pca_5five_sub5_100bp_Contri[1], "%)", sep="" )
label2_5five_sub5_100bp <-   paste( "PC2 ",  "(", myres.pca_5five_sub5_100bp_Contri[2], "%)", sep="" )
label3_5five_sub5_100bp <-   paste( "PC3 ",  "(", myres.pca_5five_sub5_100bp_Contri[3], "%)", sep="" )
label1_5five_sub5_100bp  
label2_5five_sub5_100bp  
label3_5five_sub5_100bp  


myLabel = myType1
dataframeB_5five_sub5_100bp  <- data.frame( as.data.frame(myres.pca_5five_sub5_100bp_matrix), myType2, myType3, myLabel   ) 
dataframeB_5five_sub5_100bp 


FigureTemp1_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="23A-PCA-PC1-PC2-100bpBin",  height1=3,  width1=4.8)


FigureTemp2_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="23B-PCA-PC1-PC2-alpha-100bpBin",   height1=3,  width1=4.8)


FigureTemp3_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="23C-PCA-PC1-PC2-smallDot-100bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="23D-PCA-PC1-PC2-big-100bpBin",   height1=3,  width1=4.8)



FigureTemp5_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="23E-PCA-PC1-PC2-text-100bpBin",   height1=3,  width1=4.8)


FigureTemp6_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label2_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="23F-PCA-PC1-PC2-text2-100bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="24A-PCA-PC1-PC3-100bpBin",   height1=3,  width1=4.8)


FigureTemp12_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="24B-PCA-PC1-PC3-alpha-100bpBin",   height1=3,  width1=4.8)


FigureTemp13_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="24C-PCA-PC1-PC3-smallDot-100bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="24D-PCA-PC1-PC3-big-100bpBin",  height1=3,  width1=4.8)



FigureTemp15_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="24E-PCA-PC1-PC3-text-100bpBin",   height1=3,  width1=4.8)


FigureTemp16_5five_sub5_100bp  <- ggplot( data = dataframeB_5five_sub5_100bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_5five_sub5_100bp) +   ylab(label3_5five_sub5_100bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_5five_sub5_100bp ,  path1=myOutDir_5five, fileName1="24F-PCA-PC1-PC3-text2-100bpBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_5five, "25A_PCA-3d-100bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_5five_sub5_100bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_5five_sub5_100bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_5five_sub5_100bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_5five_sub5_100bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_5five_sub5_100bp[,1:3], pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_5five_sub5_100bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_5five, "25B_PCA-3d-100bpBin-withLabels.pdf",  sep="/") )

s3d1_5five_sub5_100bp <- scatterplot3d( 
  dataframeB_5five_sub5_100bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_5five_sub5_100bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_5five_sub5_100bp$xyz.convert(dataframeB_5five_sub5_100bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_5five_sub5_100bp <- scatterplot3d( 
  dataframeB_5five_sub5_100bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_5five_sub5_100bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_5five_sub5_100bp$xyz.convert(dataframeB_5five_sub5_100bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_5five_sub5_100bp <- scatterplot3d( 
  dataframeB_5five_sub5_100bp[,1:3], pch =myType2_shape, 
  xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_5five_sub5_100bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_5five_sub5_100bp$xyz.convert(dataframeB_5five_sub5_100bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_5five, "26A_PCA-3d-100bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_5five, "26B_PCA-3d-100bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_5five, "27_PCA-3d-100bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )



##
scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )



##
scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )


scatter3D(x = dataframeB_5five_sub5_100bp[,1], y = dataframeB_5five_sub5_100bp[,2], z = dataframeB_5five_sub5_100bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_5five_sub5_100bp,  ylab = label2_5five_sub5_100bp,  zlab = label3_5five_sub5_100bp  )

dev.off()






#################



km.res1_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=2, iter.max = 20 )
km.res2_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=3, iter.max = 20  )
km.res3_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=4, iter.max = 20  )
km.res4_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=5, iter.max = 20 )
km.res5_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=6, iter.max = 20 )
km.res6_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=7, iter.max = 20 )
km.res7_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=8, iter.max = 20 )
km.res8_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=9, iter.max = 20 )
km.res9_5five_sub5_100bp  <- kmeans( t(mat_5five_sub5_100bp), centers=10, iter.max = 20 )
km.res10_5five_sub5_100bp <- kmeans( t(mat_5five_sub5_100bp), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_5five, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_5five_sub5_100bp$cluster
print("#################### 3 classes:")
km.res2_5five_sub5_100bp$cluster
print("#################### 4 classes:")
km.res3_5five_sub5_100bp$cluster
print("#################### 5 classes:")
km.res4_5five_sub5_100bp$cluster
print("#################### 6 classes:")
km.res5_5five_sub5_100bp$cluster
print("#################### 7 classes:")
km.res6_5five_sub5_100bp$cluster
print("#################### 8 classes:")
km.res7_5five_sub5_100bp$cluster
print("#################### 9 classes:")
km.res8_5five_sub5_100bp$cluster
print("#################### 10 classes:")
km.res9_5five_sub5_100bp$cluster
print("#################### 11 classes:")
km.res10_5five_sub5_100bp$cluster
sink()



pdf( file = paste(myOutDir_5five, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_5five_sub5_100bp), col = km.res1_5five_sub5_100bp$cluster )

plot( t(mat_5five_sub5_100bp), col = km.res1_5five_sub5_100bp$cluster)
points(km.res1_5five_sub5_100bp$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_5five_sub5_100bp), col = km.res1_5five_sub5_100bp$cluster)
points(km.res1_5five_sub5_100bp$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_5five_sub5_100bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 






pdf( file = paste(myOutDir_5five, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_5five_sub5_100bp), col = km.res2_5five_sub5_100bp$cluster )

plot( t(mat_5five_sub5_100bp), col = km.res2_5five_sub5_100bp$cluster)
points(km.res2_5five_sub5_100bp$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_5five_sub5_100bp), col = km.res2_5five_sub5_100bp$cluster)
points(km.res2_5five_sub5_100bp$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_5five_sub5_100bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_5five, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_5five_sub5_100bp), col = km.res3_5five_sub5_100bp$cluster )

plot( t(mat_5five_sub5_100bp), col = km.res3_5five_sub5_100bp$cluster)
points(km.res3_5five_sub5_100bp$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_5five_sub5_100bp), col = km.res3_5five_sub5_100bp$cluster)
points(km.res3_5five_sub5_100bp$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_5five_sub5_100bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 







pdf( file = paste(myOutDir_5five, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_5five_sub5_100bp), col = km.res4_5five_sub5_100bp$cluster )

plot( t(mat_5five_sub5_100bp), col = km.res4_5five_sub5_100bp$cluster)
points(km.res4_5five_sub5_100bp$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_5five_sub5_100bp), col = km.res4_5five_sub5_100bp$cluster)
points(km.res4_5five_sub5_100bp$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_5five_sub5_100bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "euclidean")
res.dist2_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "maximum")
res.dist3_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "manhattan")
res.dist4_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "canberra")
res.dist5_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "binary")
res.dist6_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "minkowski")
res.dist7_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "pearson")
res.dist8_5five_sub5_100bp <- get_dist( t(mat_5five_sub5_100bp) ,   method = "spearman")




pdf( file = paste(myOutDir_5five, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_5five_sub5_100bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_5five, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_5five_sub5_100bp ) 
print("################### maximum: ")
print(res.dist2_5five_sub5_100bp ) 
print("################### manhattan: ")
print(res.dist3_5five_sub5_100bp ) 
print("################### canberra: ")
print(res.dist4_5five_sub5_100bp ) 
print("################### binary: ")
print(res.dist5_5five_sub5_100bp ) 
print("################### minkowski: ")
print(res.dist6_5five_sub5_100bp ) 
print("################### pearson: ")
print(res.dist7_5five_sub5_100bp ) 
print("################### spearman: ")
print(res.dist8_5five_sub5_100bp ) 
sink() 


# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "ward.D" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "ward.D" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "ward.D" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "ward.D" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "ward.D" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "ward.D" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "ward.D" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "ward.D" )  

pdf( file = paste(myOutDir_5five, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "ward.D2" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "ward.D2" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "ward.D2" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "ward.D2" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "ward.D2" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "ward.D2" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "ward.D2" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "ward.D2" )  

pdf( file = paste(myOutDir_5five, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "single" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "single" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "single" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "single" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "single" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "single" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "single" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "single" )  

pdf( file = paste(myOutDir_5five, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "complete" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "complete" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "complete" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "complete" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "complete" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "complete" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "complete" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "complete" )  

pdf( file = paste(myOutDir_5five, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "average" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "average" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "average" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "average" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "average" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "average" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "average" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "average" )  

pdf( file = paste(myOutDir_5five, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "mcquitty" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "mcquitty" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "mcquitty" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "mcquitty" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "mcquitty" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "mcquitty" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "mcquitty" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "mcquitty" )  

pdf( file = paste(myOutDir_5five, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "median" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "median" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "median" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "median" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "median" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "median" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "median" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "median" )  

pdf( file = paste(myOutDir_5five, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_5five_sub5_100bp <- hclust(res.dist1_5five_sub5_100bp, method = "centroid" )   
res.hc2_5five_sub5_100bp <- hclust(res.dist2_5five_sub5_100bp, method = "centroid" )   
res.hc3_5five_sub5_100bp <- hclust(res.dist3_5five_sub5_100bp, method = "centroid" )   
res.hc4_5five_sub5_100bp <- hclust(res.dist4_5five_sub5_100bp, method = "centroid" )   
res.hc5_5five_sub5_100bp <- hclust(res.dist5_5five_sub5_100bp, method = "centroid" )   
res.hc6_5five_sub5_100bp <- hclust(res.dist6_5five_sub5_100bp, method = "centroid" )   
res.hc7_5five_sub5_100bp <- hclust(res.dist7_5five_sub5_100bp, method = "centroid" )   
res.hc8_5five_sub5_100bp <- hclust(res.dist8_5five_sub5_100bp, method = "centroid" )  

pdf( file = paste(myOutDir_5five, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_5five_sub5_100bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_5five, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp )
fviz_dend(res.hc2_5five_sub5_100bp )
fviz_dend(res.hc3_5five_sub5_100bp )
fviz_dend(res.hc4_5five_sub5_100bp )
fviz_dend(res.hc5_5five_sub5_100bp )
fviz_dend(res.hc6_5five_sub5_100bp )
fviz_dend(res.hc7_5five_sub5_100bp )
fviz_dend(res.hc8_5five_sub5_100bp )
dev.off() 


pdf( file = paste(myOutDir_5five, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc2_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc3_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc4_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc5_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc6_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc7_5five_sub5_100bp, label_cols = myType3_color )
fviz_dend(res.hc8_5five_sub5_100bp, label_cols = myType3_color )
dev.off() 









##################
myOutDir_5five_diffMe   <- paste(myOutDir_5five, "diffMe_DMC_DMR",  sep="/")
if( ! file.exists(myOutDir_5five_diffMe) ) { dir.create(myOutDir_5five_diffMe, recursive = TRUE) }


myDiff_5five_sub5_100bp = calculateDiffMeth(meth_5five_sub5_100bp,  num.cores=16)
dim(myDiff_5five_sub5_100bp)
names(myDiff_5five_sub5_100bp)
head(myDiff_5five_sub5_100bp)

myQvalue_5five_sub5_100bp = myDiff_5five_sub5_100bp$qvalue
myMethDi_5five_sub5_100bp = myDiff_5five_sub5_100bp$meth.diff



pdf(paste(myOutDir_5five_diffMe,  "1A_qvalue_distribution.pdf",  sep="/"))
hist( myQvalue_5five_sub5_100bp,  nclass=100, xlim=c(0, 1), freq=FALSE)
hist( myQvalue_5five_sub5_100bp[myQvalue_5five_sub5_100bp<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE)
hist( myQvalue_5five_sub5_100bp,  nclass=100, xlim=c(0, 1), freq=TRUE)
hist( myQvalue_5five_sub5_100bp[myQvalue_5five_sub5_100bp<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE)
dev.off() 


pdf(paste(myOutDir_5five_diffMe,  "1B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_5five_sub5_100bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_5five_sub5_100bp, nclass=100, xlim=c(0, 50), freq=FALSE)
hist(myMethDi_5five_sub5_100bp[myMethDi_5five_sub5_100bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_5five_sub5_100bp[myMethDi_5five_sub5_100bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_5five_sub5_100bp[myMethDi_5five_sub5_100bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
hist(myMethDi_5five_sub5_100bp, nclass=100, xlim=c(0, 100), freq=TRUE)
hist(myMethDi_5five_sub5_100bp, nclass=100, xlim=c(0, 50), freq=TRUE)
hist(myMethDi_5five_sub5_100bp[myMethDi_5five_sub5_100bp>=20], nclass=81, xlim=c(20, 100),  freq=TRUE)
hist(myMethDi_5five_sub5_100bp[myMethDi_5five_sub5_100bp>=10], nclass=91, xlim=c(10, 100),  freq=TRUE)
hist(myMethDi_5five_sub5_100bp[myMethDi_5five_sub5_100bp<=10], nclass=100, xlim=c(0, 10),   freq=TRUE)
dev.off() 




myColor1_5five_sub5_100bp <- rep( "no",   times= length(myQvalue_5five_sub5_100bp) )
myColor1_5five_sub5_100bp[ (abs(myMethDi_5five_sub5_100bp)>10) & (myQvalue_5five_sub5_100bp<0.001) ]  <- "yes"
length( myColor1_5five_sub5_100bp[myColor1_5five_sub5_100bp=="yes"] )



DataFrame2_5five_sub5_100bp <- data.frame(myx1 =  myMethDi_5five_sub5_100bp,   
                                          myy1 =  -log10(myQvalue_5five_sub5_100bp),  
                                          mycolor1 = myColor1_5five_sub5_100bp )

ggplot(DataFrame2_5five_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_5five_diffMe,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_5five_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_5five_diffMe,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_5five_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_5five_diffMe,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )




ggplot(DataFrame2_5five_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_5five_diffMe,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_5five_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_5five_diffMe,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_5five_sub5_100bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_5five_diffMe,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )









myDiff25p.hypo_5five_sub5_100bp  = getMethylDiff(myDiff_5five_sub5_100bp, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_5five_sub5_100bp = getMethylDiff(myDiff_5five_sub5_100bp, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
myDiff25p_5five_sub5_100bp       = getMethylDiff(myDiff_5five_sub5_100bp, difference=10, qvalue=0.05)
myDiffTemp_5five_sub5_100bp      = getMethylDiff(myDiff_5five_sub5_100bp, difference=0,  qvalue=0.05)

write.table(myDiff_5five_sub5_100bp , file = paste(myOutDir_5five_diffMe,"2A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_5five_sub5_100bp , file = paste(myOutDir_5five_diffMe,"2B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_5five_sub5_100bp , file = paste(myOutDir_5five_diffMe,"2C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_5five_sub5_100bp , file = paste(myOutDir_5five_diffMe,"2D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_5five_sub5_100bp , file = paste(myOutDir_5five_diffMe,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_5five_diffMe, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_5five_sub5_100bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_5five_diffMe, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_5five_sub5_100bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
diffGeneAnn_5five_sub5_100bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_5five_sub5_100bp,  "GRanges"),  gene.obj)

sink( file=paste(myOutDir_5five_diffMe, "4A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffGeneAnn_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_5five_sub5_100bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_5five_sub5_100bp,"GRanges"),
                                                            cpg.obj$CpGi,  cpg.obj$shores,
                                                            feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "4C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffCpGann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_5five_sub5_100bp,"GRanges"),
                                                               myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                               feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "4E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffrepeatann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_5five_sub5_100bp, "GRanges"),
                                                                   feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "5A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_5five_sub5_100bp, "GRanges"),
                                                                   feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "5C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_5five_sub5_100bp, "GRanges"),
                                                                   feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "5E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_5five_sub5_100bp, "GRanges"),
                                                                   feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "5G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_5five_sub5_100bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_5five_sub5_100bp, "GRanges"),
                                                                   feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                   feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "5I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_5five_sub5_100bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_5five_sub5_100bp_hypo)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_5five_sub5_100bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_5five_sub5_100bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_5five_sub5_100bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_5five_diffMe, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffGeneAnn_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_5five_sub5_100bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_5five_sub5_100bp,"GRanges"),
                                                             cpg.obj$CpGi,  cpg.obj$shores,
                                                             feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffCpGann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_5five_sub5_100bp,"GRanges"),
                                                                myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                                feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffrepeatann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_5five_sub5_100bp, "GRanges"),
                                                                    feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_5five_sub5_100bp, "GRanges"),
                                                                    feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_5five_sub5_100bp, "GRanges"),
                                                                    feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_5five_sub5_100bp, "GRanges"),
                                                                    feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_5five_sub5_100bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_5five_sub5_100bp, "GRanges"),
                                                                    feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                    feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_5five_diffMe, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_5five_sub5_100bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_5five_sub5_100bp_hyper)
sink()

pdf( file=paste(myOutDir_5five_diffMe, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_5five_sub5_100bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##################################################################################################################
##################################################################################################################









































######################################################################################################################################################
######################################################################################################################################################
myOutDir_6six = paste(myOutDir, "/6-1bpBin",  sep="");
if( ! file.exists(myOutDir_6six) ) { dir.create(myOutDir_6six, recursive = TRUE) }


#####
meth_6six_sub6_1bp  = unite( myobj, destrand=FALSE   )   ## 100% overlap
mat_6six_sub6_1bp   = percMethylation( meth_6six_sub6_1bp )


sink( file=paste(myOutDir_6six , "0-dimensions-1bpBin.txt", sep="/")  )
print( tiles_6six_sub6_1bp )
print("#########dimensions:")
print( dim(meth_6six_sub6_1bp)  )   
print( dim(mat_6six_sub6_1bp)   )
sink()




write.table(meth_6six_sub6_1bp , 
            file = paste(myOutDir_6six,   "1A-meth-1bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_6six_sub6_1bp , 
            file = paste(myOutDir_6six,   "1B-mat-1bpBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





pdf( file=paste(myOutDir_6six, "2A-MethylationStats-1bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getMethylationStats(tiles_6six_sub6_1bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()

sink( file=paste(myOutDir_6six, "2B-MethylationStats-1bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getMethylationStats( tiles_6six_sub6_1bp[[i]] )  )
}
sink()




pdf( file=paste(myOutDir_6six, "3A-CoverageStats-1bpBin.pdf", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  getCoverageStats(tiles_6six_sub6_1bp[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()
sink( file=paste(myOutDir_6six, "3B-CoverageStats-1bpBin.txt", sep="/")  )
for( i in c(1:length(myFileLists)) ) {
  print("##############################################")
  print( paste(mySampleID[i], "-",  myTreatment[i], ":", sep="") )
  print( getCoverageStats( tiles_6six_sub6_1bp[[i]] )  )
}
sink()





MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4A-cluster-kept100percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0 )

MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4B-cluster-kept80percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.2 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4C-cluster-kept60percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.4 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4D-cluster-kept50percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.5 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4E-cluster-kept40percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.6 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4F-cluster-kept30percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.7 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4G-cluster-kept20percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.8 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4H-cluster-kept10percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.9 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4I-cluster-kept5percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.95 )


MyCluster_1(  mymeth1 = meth_6six_sub6_1bp , path1 = myOutDir_6six, 
              file1 = "4J-cluster-kept1percent-1bpBin.pdf",   width1 = 8, height1 = 5, sdThres1 = 0.99 )





pdf( file=paste(myOutDir_6six, "5-PCA-1bpBin.pdf", sep="/") , width=8, height=5  )
PCASamples(meth_6six_sub6_1bp , screeplot=TRUE)
PCASamples(meth_6six_sub6_1bp )
dev.off()

#####################





PCA_6six_sub6_1bp <- prcomp( t(mat_6six_sub6_1bp)  )
names(PCA_6six_sub6_1bp)

sink( file = paste(myOutDir_6six,  "6A-PCA-1bpBin.txt",  sep="/") )
print(PCA_6six_sub6_1bp)
sink()

sink( file = paste(myOutDir_6six,  "6B-PCA-summary-1bpBin.txt",  sep="/") )
summary(PCA_6six_sub6_1bp)
sink()


sink( file = paste(myOutDir_6six,  "6C-PCA-all-1bpBin.txt",  sep="/") )
print("####################### myPCA2$sdev #########################")
print(PCA_6six_sub6_1bp$sdev)
print("####################### myPCA2$rotation #########################")
print(PCA_6six_sub6_1bp$rotation)
print("####################### myPCA2$center #########################")
print(PCA_6six_sub6_1bp$center)
print("####################### myPCA2$scale #########################")
print(PCA_6six_sub6_1bp$scale)
print("####################### myPCA2$x #########################")
print(PCA_6six_sub6_1bp$x)
sink()


pdf( file=paste(myOutDir_6six,   "6D-PCA-info-1bpBin.pdf", sep="/")  )
plot(PCA_6six_sub6_1bp, type="lines")
fviz_eig(PCA_6six_sub6_1bp)
dev.off() 




my_fviz_pca_ind1_6six_sub6_1bp <- fviz_pca_ind(PCA_6six_sub6_1bp,
                                               col.ind = "cos2", # Color by the quality of representation
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_6six_sub6_1bp <- fviz_pca_ind(PCA_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)
my_fviz_pca_ind3_6six_sub6_1bp <- fviz_pca_ind(PCA_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ) , # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind4_6six_sub6_1bp <- fviz_pca_ind(PCA_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               #legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)


svg(file=paste(myOutDir_6six, "7A-PCA-2D-1-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "7B-PCA-2D-2-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "7C-PCA-2D-3-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "7D-PCA-2D-4-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_6six_sub6_1bp)
dev.off() 








#############################


PCA_6six_sub6_1bp_matrix <- PCA_6six_sub6_1bp$x
dim( PCA_6six_sub6_1bp_matrix )

PCA_6six_sub6_1bp_Contri  <- (PCA_6six_sub6_1bp$sdev)^2
PCA_6six_sub6_1bp_Contri  <- PCA_6six_sub6_1bp_Contri/sum(PCA_6six_sub6_1bp_Contri)
PCA_6six_sub6_1bp_Contri  <- PCA_6six_sub6_1bp_Contri * 100
PCA_6six_sub6_1bp_Contri  <- round(PCA_6six_sub6_1bp_Contri, 2)

label1_6six_sub6_1bp <-   paste( "PC1 ",  "(", PCA_6six_sub6_1bp_Contri[1], "%)", sep="" )
label2_6six_sub6_1bp <-   paste( "PC2 ",  "(", PCA_6six_sub6_1bp_Contri[2], "%)", sep="" )
label3_6six_sub6_1bp <-   paste( "PC3 ",  "(", PCA_6six_sub6_1bp_Contri[3], "%)", sep="" )
label1_6six_sub6_1bp  
label2_6six_sub6_1bp  
label3_6six_sub6_1bp  


myLabel = myType1
dataframeA_6six_sub6_1bp  <- data.frame( as.data.frame(PCA_6six_sub6_1bp_matrix), myType2, myType3, myLabel   ) 
dataframeA_6six_sub6_1bp 


FigureTemp1_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="8A-PCA-PC1-PC2-1bpBin",  height1=3,  width1=4.8)


FigureTemp2_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="8B-PCA-PC1-PC2-alpha-1bpBin",   height1=3,  width1=4.8)


FigureTemp3_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="8C-PCA-PC1-PC2-smallDot-1bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="8D-PCA-PC1-PC2-big-1bpBin",   height1=3,  width1=4.8)



FigureTemp5_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="8E-PCA-PC1-PC2-text-1bpBin",   height1=3,  width1=4.8)


FigureTemp6_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="8F-PCA-PC1-PC2-text2-1bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="9A-PCA-PC1-PC3-1bpBin",   height1=3,  width1=4.8)


FigureTemp12_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="9B-PCA-PC1-PC3-alpha-1bpBin",   height1=3,  width1=4.8)


FigureTemp13_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="9C-PCA-PC1-PC3-smallDot-1bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="9D-PCA-PC1-PC3-big-1bpBin",  height1=3,  width1=4.8)



FigureTemp15_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="9E-PCA-PC1-PC3-text-1bpBin",   height1=3,  width1=4.8)


FigureTemp16_6six_sub6_1bp  <- ggplot( data = dataframeA_6six_sub6_1bp, aes(x = PC1, y = PC3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="9F-PCA-PC1-PC3-text2-1bpBin",   height1=3,  width1=4.8)
##################





## 3D plot: plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_6six, "10A_PCA-3d-1bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeA_6six_sub6_1bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_6six_sub6_1bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_6six_sub6_1bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_6six_sub6_1bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeA_6six_sub6_1bp[,1:3], pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_6six_sub6_1bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_6six, "10B_PCA-3d-1bpBin-withLabels.pdf",  sep="/") )

s3d1_6six_sub6_1bp <- scatterplot3d( 
  dataframeA_6six_sub6_1bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_6six_sub6_1bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_6six_sub6_1bp$xyz.convert(dataframeA_6six_sub6_1bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_6six_sub6_1bp <- scatterplot3d( 
  dataframeA_6six_sub6_1bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_6six_sub6_1bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_6six_sub6_1bp$xyz.convert(dataframeA_6six_sub6_1bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_6six_sub6_1bp <- scatterplot3d( 
  dataframeA_6six_sub6_1bp[,1:3], pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeA_6six_sub6_1bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_6six_sub6_1bp$xyz.convert(dataframeA_6six_sub6_1bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_6six, "11A_PCA-3d-1bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_6six, "11B_PCA-3d-1bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")





pdf( file = paste(myOutDir_6six, "12_PCA-3d-1bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )



##
scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )



##
scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeA_6six_sub6_1bp[,1], y = dataframeA_6six_sub6_1bp[,2], z = dataframeA_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )

dev.off()







##################
# get_eigenvalue(res.pca): Extract the eigenvalues/variances of principal components
# fviz_eig(res.pca): Visualize the eigenvalues
# get_pca_ind(res.pca), get_pca_var(res.pca): Extract the results for individuals and variables, respectively.
# fviz_pca_ind(res.pca), fviz_pca_var(res.pca): Visualize the results individuals and variables, respectively.
# fviz_pca_biplot(res.pca): Make a biplot of individuals and variables.


myres.pca_6six_sub6_1bp  <- PCA(t(mat_6six_sub6_1bp),   graph = FALSE)

sink( file = paste(myOutDir_6six, "20_PCA-info-1bpBin-byPCA.txt",  sep="/") )
print( myres.pca_6six_sub6_1bp )
print( "#################################" )
print( summary(myres.pca_6six_sub6_1bp) )
print( "#################################"  )
myeig.val_6six_sub6_1bp <- get_eigenvalue( myres.pca_6six_sub6_1bp )
myeig.val_6six_sub6_1bp
sink() 


pdf( file = paste(myOutDir_6six, "21_PCA-screePlot-1bpBin-byPCA.pdf",  sep="/") )
fviz_eig(myres.pca_6six_sub6_1bp, addlabels = TRUE )
fviz_screeplot(X=myres.pca_6six_sub6_1bp, choice = "variance", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_6six_sub6_1bp, choice = "eigenvalue", geom = "line",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_6six_sub6_1bp, choice = "variance", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
fviz_screeplot(X=myres.pca_6six_sub6_1bp, choice = "eigenvalue", geom = "bar",
               barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
               ncp = 10, addlabels = TRUE)
dev.off()





my_fviz_pca_ind1_6six_sub6_1bp <- fviz_pca_ind(myres.pca_6six_sub6_1bp,
                                               col.ind = "cos2", # Color by the quality of representation
                                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                               repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2_6six_sub6_1bp <- fviz_pca_ind(myres.pca_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)
my_fviz_pca_ind3_6six_sub6_1bp <- fviz_pca_ind(myres.pca_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ) , # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               addEllipses = TRUE, # Concentration ellipses
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind4_6six_sub6_1bp <- fviz_pca_ind(myres.pca_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               #legend.title = "Groups",
                                               repel = TRUE, 
                                               label = "none", 
                                               alpha.ind = 1
)
my_fviz_pca_ind5_6six_sub6_1bp <- fviz_pca_ind(myres.pca_6six_sub6_1bp,
                                               col.ind =  as.factor( myTreatment ), # color by groups
                                               #palette = c("#00AFBB",  "#FC4E07"),
                                               ellipse.type = "confidence",
                                               legend.title = "Groups",
                                               repel = TRUE
)

svg(file=paste(myOutDir_6six, "22A-PCA-2D-1-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind1_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "22B-PCA-2D-2-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind2_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "22C-PCA-2D-3-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind3_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "22D-PCA-2D-4-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind4_6six_sub6_1bp)
dev.off() 

svg(file=paste(myOutDir_6six, "22E-PCA-2D-4-1bpBin.svg", sep="/") )
print(my_fviz_pca_ind5_6six_sub6_1bp)
dev.off() 



#############################


myres.pca_6six_sub6_1bp_matrix <- myres.pca_6six_sub6_1bp$ind$coord
dim( myres.pca_6six_sub6_1bp_matrix )

myres.pca_6six_sub6_1bp_Contri  <- (myres.pca_6six_sub6_1bp$eig)[,2] 
myres.pca_6six_sub6_1bp_Contri  <- myres.pca_6six_sub6_1bp_Contri/sum(myres.pca_6six_sub6_1bp_Contri)
myres.pca_6six_sub6_1bp_Contri  <- myres.pca_6six_sub6_1bp_Contri * 100
myres.pca_6six_sub6_1bp_Contri  <- round(myres.pca_6six_sub6_1bp_Contri, 2)

label1_6six_sub6_1bp <-   paste( "PC1 ",  "(", myres.pca_6six_sub6_1bp_Contri[1], "%)", sep="" )
label2_6six_sub6_1bp <-   paste( "PC2 ",  "(", myres.pca_6six_sub6_1bp_Contri[2], "%)", sep="" )
label3_6six_sub6_1bp <-   paste( "PC3 ",  "(", myres.pca_6six_sub6_1bp_Contri[3], "%)", sep="" )
label1_6six_sub6_1bp  
label2_6six_sub6_1bp  
label3_6six_sub6_1bp  


myLabel = myType1
dataframeB_6six_sub6_1bp  <- data.frame( as.data.frame(myres.pca_6six_sub6_1bp_matrix), myType2, myType3, myLabel   ) 
dataframeB_6six_sub6_1bp 


FigureTemp1_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="23A-PCA-PC1-PC2-1bpBin",  height1=3,  width1=4.8)


FigureTemp2_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="23B-PCA-PC1-PC2-alpha-1bpBin",   height1=3,  width1=4.8)


FigureTemp3_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="23C-PCA-PC1-PC2-smallDot-1bpBin",   height1=3,  width1=4.8)


FigureTemp4 <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="23D-PCA-PC1-PC2-big-1bpBin",   height1=3,  width1=4.8)



FigureTemp5_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp5_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="23E-PCA-PC1-PC2-text-1bpBin",   height1=3,  width1=4.8)


FigureTemp6_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.2,  shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label2_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp6_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="23F-PCA-PC1-PC2-text2-1bpBin",   height1=3,  width1=4.8)






## plot for PC1 and PC3
FigureTemp11_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.3,  shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp11_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="24A-PCA-PC1-PC3-1bpBin",   height1=3,  width1=4.8)


FigureTemp12_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=0.5  )+ xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp12_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="24B-PCA-PC1-PC3-alpha-1bpBin",   height1=3,  width1=4.8)


FigureTemp13_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=2, alpha=1  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp13_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="24C-PCA-PC1-PC3-smallDot-1bpBin",   height1=3,  width1=4.8)


FigureTemp14 <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=6, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +     
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp14_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="24D-PCA-PC1-PC3-big-1bpBin",  height1=3,  width1=4.8)



FigureTemp15_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=2, alpha=0.7  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp15_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="24E-PCA-PC1-PC3-text-1bpBin",   height1=3,  width1=4.8)


FigureTemp16_6six_sub6_1bp  <- ggplot( data = dataframeB_6six_sub6_1bp, aes(x = Dim.1, y = Dim.3, shape=as.factor(myType2), color=as.factor(myType3), label=myLabel )) + 
  geom_point(size=5, alpha=0.5  ) + xlab(label1_6six_sub6_1bp) +   ylab(label3_6six_sub6_1bp) +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp16_6six_sub6_1bp ,  path1=myOutDir_6six, fileName1="24F-PCA-PC1-PC3-text2-1bpBin",   height1=3,  width1=4.8)
##################





## 3D plot:  scatter3D(),  plot3d(), scatter3d(), scatterplot3d(),  ggplot2 and Plotly
##################
pdf( file = paste(myOutDir_6six, "25A_PCA-3d-1bpBin.pdf",  sep="/") )

scatterplot3d( 
  dataframeB_6six_sub6_1bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_6six_sub6_1bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_6six_sub6_1bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_6six_sub6_1bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)


scatterplot3d( 
  dataframeB_6six_sub6_1bp[,1:3], pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_6six_sub6_1bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)


dev.off()





pdf( file = paste(myOutDir_6six, "25B_PCA-3d-1bpBin-withLabels.pdf",  sep="/") )

s3d1_6six_sub6_1bp <- scatterplot3d( 
  dataframeB_6six_sub6_1bp[,1:3], color=myType3_color,   pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_6six_sub6_1bp$myType3),
       col = myType3_color2,  pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d1_6six_sub6_1bp$xyz.convert(dataframeB_6six_sub6_1bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d2_6six_sub6_1bp <- scatterplot3d( 
  dataframeB_6six_sub6_1bp[,1:3], color=myType3_color, pch=19,   
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_6six_sub6_1bp$myType3),
       col = myType3_color2, pch=19,  inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d2_6six_sub6_1bp$xyz.convert(dataframeB_6six_sub6_1bp[, 1:3]), labels = myLabel,  cex= 1 )


s3d3_6six_sub6_1bp <- scatterplot3d( 
  dataframeB_6six_sub6_1bp[,1:3], pch =myType2_shape, 
  xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
  cex.symbols = 2,  cex.axis=1, cex.lab=1.5
)
legend("top", legend = levels(dataframeB_6six_sub6_1bp$myType2),
       pch = myType2_shape2,   inset = -0.1, xpd = TRUE, horiz = TRUE)
text(s3d3_6six_sub6_1bp$xyz.convert(dataframeB_6six_sub6_1bp[, 1:3]), labels = myLabel,  cex= 1 )

dev.off()
##################





scatter3d(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_6six, "26A_PCA-3d-1bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")

scatter3d(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          groups = as.factor(myType3),   grid = FALSE,  surface = FALSE,  
          ellipsoid = FALSE,  surface.col = myType3_color2,   axis.scales = FALSE, 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp, 
          axis.col = c("black", "black", "black")
)
rgl.postscript( paste(myOutDir_6six, "26B_PCA-3d-1bpBin-byScatter3d.pdf",  sep="/"),  fmt="pdf")






pdf( file = paste(myOutDir_6six, "27_PCA-3d-1bpBin-by-scatter3D.pdf",  sep="/") )
##

##
scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30,  
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20,  
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20,   
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10,  
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30,   
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20,   
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )



##
scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )



##
scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=0, phi=30, bty = "g" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=20, phi=20, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=30, phi=20, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=10, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=45, phi=30, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )


scatter3D(x = dataframeB_6six_sub6_1bp[,1], y = dataframeB_6six_sub6_1bp[,2], z = dataframeB_6six_sub6_1bp[,3], 
          colvar = NULL, col = myType3_color,  pch = myType2_shape,   cex = 2, theta=60, phi=20, bty = "b2" , 
          xlab = label1_6six_sub6_1bp,  ylab = label2_6six_sub6_1bp,  zlab = label3_6six_sub6_1bp  )

dev.off()






#################



km.res1_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=2, iter.max = 20 )
km.res2_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=3, iter.max = 20  )
km.res3_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=4, iter.max = 20  )
km.res4_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=5, iter.max = 20 )
km.res5_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=6, iter.max = 20 )
km.res6_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=7, iter.max = 20 )
km.res7_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=8, iter.max = 20 )
km.res8_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=9, iter.max = 20 )
km.res9_6six_sub6_1bp  <- kmeans( t(mat_6six_sub6_1bp), centers=10, iter.max = 20 )
km.res10_6six_sub6_1bp <- kmeans( t(mat_6six_sub6_1bp), centers=11, iter.max = 20 )


sink( file = paste(myOutDir_6six, "50_kmeans-cluster.txt",  sep="/") )
print("#################### 2 classes:")
km.res1_6six_sub6_1bp$cluster
print("#################### 3 classes:")
km.res2_6six_sub6_1bp$cluster
print("#################### 4 classes:")
km.res3_6six_sub6_1bp$cluster
print("#################### 5 classes:")
km.res4_6six_sub6_1bp$cluster
print("#################### 6 classes:")
km.res5_6six_sub6_1bp$cluster
print("#################### 7 classes:")
km.res6_6six_sub6_1bp$cluster
print("#################### 8 classes:")
km.res7_6six_sub6_1bp$cluster
print("#################### 9 classes:")
km.res8_6six_sub6_1bp$cluster
print("#################### 10 classes:")
km.res9_6six_sub6_1bp$cluster
print("#################### 11 classes:")
km.res10_6six_sub6_1bp$cluster
sink()



pdf( file = paste(myOutDir_6six, "51A_kmeans-2classes.pdf",  sep="/") )

plot( t(mat_6six_sub6_1bp), col = km.res1_6six_sub6_1bp$cluster )

plot( t(mat_6six_sub6_1bp), col = km.res1_6six_sub6_1bp$cluster)
points(km.res1_6six_sub6_1bp$centers, col = 1:2, pch = 8, cex = 2) 

plot( t(mat_6six_sub6_1bp), col = km.res1_6six_sub6_1bp$cluster)
points(km.res1_6six_sub6_1bp$centers, col = 1:2, pch = 8, cex = 2)
text(t(mat_6six_sub6_1bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 






pdf( file = paste(myOutDir_6six, "51B_kmeans-3classes.pdf",  sep="/") )

plot( t(mat_6six_sub6_1bp), col = km.res2_6six_sub6_1bp$cluster )

plot( t(mat_6six_sub6_1bp), col = km.res2_6six_sub6_1bp$cluster)
points(km.res2_6six_sub6_1bp$centers, col = 1:3, pch = 8, cex = 2) 

plot( t(mat_6six_sub6_1bp), col = km.res2_6six_sub6_1bp$cluster)
points(km.res2_6six_sub6_1bp$centers, col = 1:3, pch = 8, cex = 2)
text(t(mat_6six_sub6_1bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 




pdf( file = paste(myOutDir_6six, "51C_kmeans-4classes.pdf",  sep="/") )

plot( t(mat_6six_sub6_1bp), col = km.res3_6six_sub6_1bp$cluster )

plot( t(mat_6six_sub6_1bp), col = km.res3_6six_sub6_1bp$cluster)
points(km.res3_6six_sub6_1bp$centers, col = 1:4, pch = 8, cex = 2) 

plot( t(mat_6six_sub6_1bp), col = km.res3_6six_sub6_1bp$cluster)
points(km.res3_6six_sub6_1bp$centers, col = 1:4, pch = 8, cex = 2)
text(t(mat_6six_sub6_1bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 







pdf( file = paste(myOutDir_6six, "51D_kmeans-5classes.pdf",  sep="/") )

plot( t(mat_6six_sub6_1bp), col = km.res4_6six_sub6_1bp$cluster )

plot( t(mat_6six_sub6_1bp), col = km.res4_6six_sub6_1bp$cluster)
points(km.res4_6six_sub6_1bp$centers, col = 1:5, pch = 8, cex = 2) 

plot( t(mat_6six_sub6_1bp), col = km.res4_6six_sub6_1bp$cluster)
points(km.res4_6six_sub6_1bp$centers, col = 1:5, pch = 8, cex = 2)
text(t(mat_6six_sub6_1bp), labels=myType1, cex= 0.7, pos=1)

dev.off() 





######################

res.dist1_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "euclidean")
res.dist2_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "maximum")
res.dist3_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "manhattan")
res.dist4_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "canberra")
res.dist5_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "binary")
res.dist6_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "minkowski")
res.dist7_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "pearson")
res.dist8_6six_sub6_1bp <- get_dist( t(mat_6six_sub6_1bp) ,   method = "spearman")




pdf( file = paste(myOutDir_6six, "100A_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist5_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist6_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist7_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist8_6six_sub6_1bp,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


sink( file = paste(myOutDir_6six, "100B_all-distance-matrix.txt",  sep="/") )
print("################### euclidean: ")
print(res.dist1_6six_sub6_1bp ) 
print("################### maximum: ")
print(res.dist2_6six_sub6_1bp ) 
print("################### manhattan: ")
print(res.dist3_6six_sub6_1bp ) 
print("################### canberra: ")
print(res.dist4_6six_sub6_1bp ) 
print("################### binary: ")
print(res.dist5_6six_sub6_1bp ) 
print("################### minkowski: ")
print(res.dist6_6six_sub6_1bp ) 
print("################### pearson: ")
print(res.dist7_6six_sub6_1bp ) 
print("################### spearman: ")
print(res.dist8_6six_sub6_1bp ) 
sink() 


# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "ward.D" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "ward.D" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "ward.D" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "ward.D" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "ward.D" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "ward.D" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "ward.D" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "ward.D" )  

pdf( file = paste(myOutDir_6six, "101A_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "101B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "101C_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 







# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "ward.D2" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "ward.D2" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "ward.D2" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "ward.D2" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "ward.D2" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "ward.D2" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "ward.D2" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "ward.D2" )  

pdf( file = paste(myOutDir_6six, "102A_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "102B_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "102C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 













# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "single" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "single" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "single" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "single" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "single" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "single" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "single" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "single" )  

pdf( file = paste(myOutDir_6six, "103A_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "103B_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "103C_hierarchical-single.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "complete" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "complete" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "complete" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "complete" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "complete" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "complete" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "complete" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "complete" )  

pdf( file = paste(myOutDir_6six, "104A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "104B_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "104C_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "average" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "average" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "average" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "average" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "average" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "average" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "average" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "average" )  

pdf( file = paste(myOutDir_6six, "105A_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "105B_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "105C_hierarchical-average.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 










# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "mcquitty" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "mcquitty" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "mcquitty" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "mcquitty" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "mcquitty" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "mcquitty" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "mcquitty" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "mcquitty" )  

pdf( file = paste(myOutDir_6six, "106A_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "106B_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "106C_hierarchical-mcquitty.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 









# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "median" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "median" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "median" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "median" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "median" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "median" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "median" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "median" )  

pdf( file = paste(myOutDir_6six, "107A_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "107B_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "107C_hierarchical-median.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 








# Compute hierarchical clustering
res.hc1_6six_sub6_1bp <- hclust(res.dist1_6six_sub6_1bp, method = "centroid" )   
res.hc2_6six_sub6_1bp <- hclust(res.dist2_6six_sub6_1bp, method = "centroid" )   
res.hc3_6six_sub6_1bp <- hclust(res.dist3_6six_sub6_1bp, method = "centroid" )   
res.hc4_6six_sub6_1bp <- hclust(res.dist4_6six_sub6_1bp, method = "centroid" )   
res.hc5_6six_sub6_1bp <- hclust(res.dist5_6six_sub6_1bp, method = "centroid" )   
res.hc6_6six_sub6_1bp <- hclust(res.dist6_6six_sub6_1bp, method = "centroid" )   
res.hc7_6six_sub6_1bp <- hclust(res.dist7_6six_sub6_1bp, method = "centroid" )   
res.hc8_6six_sub6_1bp <- hclust(res.dist8_6six_sub6_1bp, method = "centroid" )  

pdf( file = paste(myOutDir_6six, "108A_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc5_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc6_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc7_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc8_6six_sub6_1bp, k = 2,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 

pdf( file = paste(myOutDir_6six, "108B_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp )
fviz_dend(res.hc2_6six_sub6_1bp )
fviz_dend(res.hc3_6six_sub6_1bp )
fviz_dend(res.hc4_6six_sub6_1bp )
fviz_dend(res.hc5_6six_sub6_1bp )
fviz_dend(res.hc6_6six_sub6_1bp )
fviz_dend(res.hc7_6six_sub6_1bp )
fviz_dend(res.hc8_6six_sub6_1bp )
dev.off() 


pdf( file = paste(myOutDir_6six, "108C_hierarchical-centroid.pdf",  sep="/") )
fviz_dend(res.hc1_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc2_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc3_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc4_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc5_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc6_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc7_6six_sub6_1bp, label_cols = myType3_color )
fviz_dend(res.hc8_6six_sub6_1bp, label_cols = myType3_color )
dev.off() 









##################
myOutDir_6six_diffMe   <- paste(myOutDir_6six, "diffMe_DMC_DMR",  sep="/")
if( ! file.exists(myOutDir_6six_diffMe) ) { dir.create(myOutDir_6six_diffMe, recursive = TRUE) }


myDiff_6six_sub6_1bp = calculateDiffMeth(meth_6six_sub6_1bp,  num.cores=16)
dim(myDiff_6six_sub6_1bp)
names(myDiff_6six_sub6_1bp)
head(myDiff_6six_sub6_1bp)

myQvalue_6six_sub6_1bp = myDiff_6six_sub6_1bp$qvalue
myMethDi_6six_sub6_1bp = myDiff_6six_sub6_1bp$meth.diff



pdf(paste(myOutDir_6six_diffMe,  "1A_qvalue_distribution.pdf",  sep="/"))
hist( myQvalue_6six_sub6_1bp,  nclass=100, xlim=c(0, 1), freq=FALSE)
hist( myQvalue_6six_sub6_1bp[myQvalue_6six_sub6_1bp<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE)
hist( myQvalue_6six_sub6_1bp,  nclass=100, xlim=c(0, 1), freq=TRUE)
hist( myQvalue_6six_sub6_1bp[myQvalue_6six_sub6_1bp<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE)
dev.off() 


pdf(paste(myOutDir_6six_diffMe,  "1B_methDiff_distribution.pdf",  sep="/"))
hist(myMethDi_6six_sub6_1bp, nclass=100, xlim=c(0, 100), freq=FALSE)
hist(myMethDi_6six_sub6_1bp, nclass=100, xlim=c(0, 50), freq=FALSE)
hist(myMethDi_6six_sub6_1bp[myMethDi_6six_sub6_1bp>=20], nclass=81, xlim=c(20, 100),  freq=FALSE)
hist(myMethDi_6six_sub6_1bp[myMethDi_6six_sub6_1bp>=10], nclass=91, xlim=c(10, 100),  freq=FALSE)
hist(myMethDi_6six_sub6_1bp[myMethDi_6six_sub6_1bp<=10], nclass=100, xlim=c(0, 10),   freq=FALSE)
hist(myMethDi_6six_sub6_1bp, nclass=100, xlim=c(0, 100), freq=TRUE)
hist(myMethDi_6six_sub6_1bp, nclass=100, xlim=c(0, 50), freq=TRUE)
hist(myMethDi_6six_sub6_1bp[myMethDi_6six_sub6_1bp>=20], nclass=81, xlim=c(20, 100),  freq=TRUE)
hist(myMethDi_6six_sub6_1bp[myMethDi_6six_sub6_1bp>=10], nclass=91, xlim=c(10, 100),  freq=TRUE)
hist(myMethDi_6six_sub6_1bp[myMethDi_6six_sub6_1bp<=10], nclass=100, xlim=c(0, 10),   freq=TRUE)
dev.off() 




myColor1_6six_sub6_1bp <- rep( "no",   times= length(myQvalue_6six_sub6_1bp) )
myColor1_6six_sub6_1bp[ (abs(myMethDi_6six_sub6_1bp)>10) & (myQvalue_6six_sub6_1bp<0.001) ]  <- "yes"
length( myColor1_6six_sub6_1bp[myColor1_6six_sub6_1bp=="yes"] )



DataFrame2_6six_sub6_1bp <- data.frame(myx1 =  myMethDi_6six_sub6_1bp,   
                                       myy1 =  -log10(myQvalue_6six_sub6_1bp),  
                                       mycolor1 = myColor1_6six_sub6_1bp )

ggplot(DataFrame2_6six_sub6_1bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_6six_diffMe,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_6six_sub6_1bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_6six_diffMe,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_6six_sub6_1bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.2 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_6six_diffMe,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )




ggplot(DataFrame2_6six_sub6_1bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  
ggsave( filename = paste(myOutDir_6six_diffMe,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_6six_sub6_1bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-40, 40)  
ggsave( filename = paste(myOutDir_6six_diffMe,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )

ggplot(DataFrame2_6six_sub6_1bp, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
  geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
  xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1(textSize1=14)  + xlim(-20, 20)  
ggsave( filename = paste(myOutDir_6six_diffMe,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )









myDiff25p.hypo_6six_sub6_1bp  = getMethylDiff(myDiff_6six_sub6_1bp, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
myDiff25p.hyper_6six_sub6_1bp = getMethylDiff(myDiff_6six_sub6_1bp, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
myDiff25p_6six_sub6_1bp       = getMethylDiff(myDiff_6six_sub6_1bp, difference=10, qvalue=0.05)
myDiffTemp_6six_sub6_1bp      = getMethylDiff(myDiff_6six_sub6_1bp, difference=0,  qvalue=0.05)

write.table(myDiff_6six_sub6_1bp , file = paste(myOutDir_6six_diffMe,"2A_diffMe-allsites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hypo_6six_sub6_1bp , file = paste(myOutDir_6six_diffMe,"2B_diffMe-hypo.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p.hyper_6six_sub6_1bp , file = paste(myOutDir_6six_diffMe,"2C_diffMe-hyper.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiff25p_6six_sub6_1bp , file = paste(myOutDir_6six_diffMe,"2D_AlldiffMesites.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

write.table(myDiffTemp_6six_sub6_1bp , file = paste(myOutDir_6six_diffMe,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")





sink( file=paste(myOutDir_6six_diffMe, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
diffMethPerChr(myDiff_6six_sub6_1bp,  plot=FALSE,  qvalue.cutoff=0.01, meth.cutoff=10)
sink()


pdf( file=paste(myOutDir_6six_diffMe, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
diffMethPerChr(myDiff_6six_sub6_1bp,  plot=TRUE,  qvalue.cutoff=0.01, meth.cutoff=10)
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
diffGeneAnn_6six_sub6_1bp_hypo = annotateWithGeneParts(as(myDiff25p.hypo_6six_sub6_1bp,  "GRanges"),  gene.obj)

sink( file=paste(myOutDir_6six_diffMe, "4A-distribution-onGenes-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffGeneAnn_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_6six_sub6_1bp_hypo,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_6six_sub6_1bp,"GRanges"),
                                                         cpg.obj$CpGi,  cpg.obj$shores,
                                                         feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "4C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffCpGann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_6six_sub6_1bp,"GRanges"),
                                                            myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                            feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "4E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffrepeatann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_6six_sub6_1bp, "GRanges"),
                                                                feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "5A-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffImprinted1Ann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_6six_sub6_1bp, "GRanges"),
                                                                feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "5C-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffImprinted2Ann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_6six_sub6_1bp, "GRanges"),
                                                                feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "5E-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffImprinted3Ann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_6six_sub6_1bp, "GRanges"),
                                                                feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "5G-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffImprinted4Ann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_6six_sub6_1bp_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_6six_sub6_1bp, "GRanges"),
                                                                feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "5I-distribution-on-hypo.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_6six_sub6_1bp_hypo,  percentage=TRUE)
print(diffImprinted5Ann_6six_sub6_1bp_hypo)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_6six_sub6_1bp_hypo, precedence=TRUE, main="differential methylation annotation")
dev.off()
####################






########################## annotation for hyper sites.
diffGeneAnn_6six_sub6_1bp_hyper = annotateWithGeneParts(as(myDiff25p.hyper_6six_sub6_1bp,"GRanges"),  gene.obj)

sink( file=paste(myOutDir_6six_diffMe, "10A-distribution-onGenes-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffGeneAnn_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffGeneAnn_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffGeneAnn_6six_sub6_1bp_hyper,precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffCpGann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_6six_sub6_1bp,"GRanges"),
                                                          cpg.obj$CpGi,  cpg.obj$shores,
                                                          feature.name="CpGi",flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "10C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffCpGann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffCpGann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffCpGann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffrepeatann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_6six_sub6_1bp,"GRanges"),
                                                             myrepeat.obj$Repeats,  myrepeat.obj$shores,
                                                             feature.name="Repeats",flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "10E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffrepeatann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffrepeatann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffrepeatann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted1Ann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_6six_sub6_1bp, "GRanges"),
                                                                 feature=imprint1.obj$ImprintedRegions,  flank=imprint1.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "11A-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted1Ann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffImprinted1Ann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted1Ann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted2Ann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_6six_sub6_1bp, "GRanges"),
                                                                 feature=imprint2.obj$ImprintedRegions,  flank=imprint2.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "11C-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted2Ann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffImprinted2Ann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted2Ann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted3Ann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_6six_sub6_1bp, "GRanges"),
                                                                 feature=imprint3.obj$ImprintedRegions,  flank=imprint3.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "11E-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted3Ann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffImprinted3Ann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted3Ann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##
diffImprinted4Ann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_6six_sub6_1bp, "GRanges"),
                                                                 feature=imprint4.obj$ImprintedRegions,  flank=imprint4.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "11G-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted4Ann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffImprinted4Ann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted4Ann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()




##
diffImprinted5Ann_6six_sub6_1bp_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_6six_sub6_1bp, "GRanges"),
                                                                 feature=imprint5.obj$ImprintedRegions,  flank=imprint5.obj$shores,
                                                                 feature.name="ImprintedRegions",  flank.name="shores")

sink( file=paste(myOutDir_6six_diffMe, "11I-distribution-on-hyper.txt", sep="/")   )
getFeatsWithTargetsStats(diffImprinted5Ann_6six_sub6_1bp_hyper,  percentage=TRUE)
print(diffImprinted5Ann_6six_sub6_1bp_hyper)
sink()

pdf( file=paste(myOutDir_6six_diffMe, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
plotTargetAnnotation(diffImprinted5Ann_6six_sub6_1bp_hyper, precedence=TRUE, main="differential methylation annotation")
dev.off()



##################################################################################################################
##################################################################################################################














