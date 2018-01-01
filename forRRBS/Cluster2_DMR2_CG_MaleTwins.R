##################################################################################################################
## All samples must be overlapped. (100% overlap)
## Suffixes of all self-defined global variables must be "_g".
##
## Example:  
## Rscript  Cluster2_DMR2_CG_MaleTwins.R     11B_splitXY/100A-rmXY/5_cov50reads     1011B_splitXY/100A-rmXY/5_cov50reads/Cluster2_DMR2_CG_MaleTwins          
         
args_g <- commandArgs(TRUE)
print("args: ")
print(args_g[1])   
print(args_g[2])     
print("#############")

inputDir_g = args_g[1];     ## the path of input files
outDir_g   = args_g[2];     ## the path of output files
# inputDir_g =  "11B_splitXY/100A-rmXY/0_cov5reads"
# outDir_g   =  "1011B_splitXY/100A-rmXY/0_cov5reads/Cluster2_DMR2_CG_MaleTwins"
print(inputDir_g)   
print(outDir_g)

if( ! file.exists(inputDir_g) ) { print("##### Error-1 #####") }
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }
##################################################################################################################





##################################################################################################################
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
library(DSS) 
library(matrixStats)
library(ggpubr)




MyTheme_1_g <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
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
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	    ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	    ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	            ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
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
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	  ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	  ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 


MySaveGgplot2_1_g <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
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


MyCluster_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  dev.off()
}


MyCluster_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=FALSE,   plot=TRUE )
  dev.off()
}


MyCluster_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995)                    
  MyCluster_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999)                    
  
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    )                     
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 )                    
  MyCluster_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.35.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.35 )                      
}


MyPCA_1_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  PCASamples(mymeth1, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE )
  dev.off()
}


MyPCA_2_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  PCASamples(mymeth1, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  PCASamples(mymeth1, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE )
  dev.off()
}


MyPCA_3_g <- function(  mymeth2 ,  path2,   file2, width2, height2 ) {
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   )                     
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 )                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995)                    
  MyPCA_1_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999)                    
  
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    )                     
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 )                    
  MyPCA_2_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.35.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.35 )                      
}



myHierarchicalClustering_1_g  <- function(  mat_3three,   path_temp1,   dataFrame_temp1  )  {
  res.dist1_3three <- get_dist( t(mat_3three) ,   method = "euclidean")
  res.dist2_3three <- get_dist( t(mat_3three) ,   method = "maximum"  )
  res.dist3_3three <- get_dist( t(mat_3three) ,   method = "manhattan")
  res.dist4_3three <- get_dist( t(mat_3three) ,   method = "canberra" )
  res.dist5_3three <- get_dist( t(mat_3three) ,   method = "binary"   )
  res.dist6_3three <- get_dist( t(mat_3three) ,   method = "minkowski")
  res.dist7_3three <- get_dist( t(mat_3three) ,   method = "pearson"  )
  res.dist8_3three <- get_dist( t(mat_3three) ,   method = "spearman" )
  
  pdf( file = paste(path_temp1, "1A_visualizing-distance-matrix.pdf",  sep="/") )
  print( fviz_dist(res.dist1_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) )
  print( fviz_dist(res.dist2_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist3_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist4_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist5_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist6_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist7_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  print( fviz_dist(res.dist8_3three,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) ) 
  dev.off() 
  
  sink( file = paste(path_temp1, "1B_all-distance-matrix.txt",  sep="/") )
  print("################### euclidean: ")
  print(res.dist1_3three ) 
  print("################### maximum: ")
  print(res.dist2_3three ) 
  print("################### manhattan: ")
  print(res.dist3_3three ) 
  print("################### canberra: ")
  print(res.dist4_3three ) 
  print("################### binary: ")
  print(res.dist5_3three ) 
  print("################### minkowski: ")
  print(res.dist6_3three ) 
  print("################### pearson: ")
  print(res.dist7_3three ) 
  print("################### spearman: ")
  print(res.dist8_3three ) 
  sink() 
  
  # Compute hierarchical clustering by "ward.D"
  res.hc1_3three <- hclust(res.dist1_3three, method = "ward.D"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "ward.D"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "ward.D"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "ward.D"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "ward.D"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "ward.D"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "ward.D"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "ward.D"  ) 
  
  pdf( file = paste(path_temp1, "2A_hierarchical-ward.D-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="ward.D, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="ward.D, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="ward.D, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="ward.D, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="ward.D, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="ward.D, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="ward.D, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="ward.D, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "2B_hierarchical-ward.D-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="ward.D, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="ward.D, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="ward.D, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="ward.D, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="ward.D, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="ward.D, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="ward.D, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="ward.D, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  # Compute hierarchical clustering by "ward.D2"
  res.hc1_3three <- hclust(res.dist1_3three, method = "ward.D2"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "ward.D2"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "ward.D2"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "ward.D2"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "ward.D2"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "ward.D2"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "ward.D2"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "ward.D2"  ) 
  
  pdf( file = paste(path_temp1, "3A_hierarchical-ward.D2-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="ward.D2, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="ward.D2, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="ward.D2, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="ward.D2, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="ward.D2, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="ward.D2, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="ward.D2, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="ward.D2, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "3B_hierarchical-ward.D2-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="ward.D2, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="ward.D2, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="ward.D2, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="ward.D2, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="ward.D2, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="ward.D2, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="ward.D2, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="ward.D2, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  # Compute hierarchical clustering by "single"
  res.hc1_3three <- hclust(res.dist1_3three, method = "single"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "single"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "single"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "single"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "single"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "single"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "single"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "single"  ) 
  
  pdf( file = paste(path_temp1, "4A_hierarchical-single-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="single, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="single, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="single, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="single, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="single, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="single, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="single, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="single, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "4B_hierarchical-single-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="single, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="single, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="single, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="single, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="single, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="single, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="single, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="single, spearman"   , repel = TRUE)  )         
  dev.off() 
  
  # Compute hierarchical clustering by "complete"
  res.hc1_3three <- hclust(res.dist1_3three, method = "complete"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "complete"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "complete"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "complete"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "complete"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "complete"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "complete"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "complete"  ) 
  
  pdf( file = paste(path_temp1, "5A_hierarchical-complete-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="complete, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="complete, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="complete, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="complete, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="complete, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="complete, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="complete, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="complete, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "5B_hierarchical-complete-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="complete, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="complete, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="complete, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="complete, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="complete, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="complete, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="complete, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="complete, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  # Compute hierarchical clustering by "average"
  res.hc1_3three <- hclust(res.dist1_3three, method = "average"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "average"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "average"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "average"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "average"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "average"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "average"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "average"  ) 
  
  pdf( file = paste(path_temp1, "6A_hierarchical-average-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="average, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="average, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="average, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="average, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="average, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="average, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="average, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="average, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "6B_hierarchical-average-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="average, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="average, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="average, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="average, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="average, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="average, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="average, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="average, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  # Compute hierarchical clustering by "mcquitty"
  res.hc1_3three <- hclust(res.dist1_3three, method = "mcquitty"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "mcquitty"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "mcquitty"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "mcquitty"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "mcquitty"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "mcquitty"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "mcquitty"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "mcquitty"  ) 
  
  pdf( file = paste(path_temp1, "7A_hierarchical-mcquitty-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="mcquitty, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="mcquitty, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="mcquitty, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="mcquitty, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="mcquitty, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="mcquitty, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="mcquitty, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="mcquitty, spearman"   )  )
  dev.off() 
  
  pdf( file = paste(path_temp1, "7B_hierarchical-mcquitty-phylogenic.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="phylogenic", main="mcquitty, euclidean"  , repel = TRUE)  )
  print(  fviz_dend(res.hc2_3three,   type="phylogenic", main="mcquitty, maximum"    , repel = TRUE)  )
  print(  fviz_dend(res.hc3_3three,   type="phylogenic", main="mcquitty, manhattan"  , repel = TRUE)  )
  print(  fviz_dend(res.hc4_3three,   type="phylogenic", main="mcquitty, canberra"   , repel = TRUE)  )
  print(  fviz_dend(res.hc5_3three,   type="phylogenic", main="mcquitty, binary"     , repel = TRUE)  )
  print(  fviz_dend(res.hc6_3three,   type="phylogenic", main="mcquitty, minkowski"  , repel = TRUE)  )
  print(  fviz_dend(res.hc7_3three,   type="phylogenic", main="mcquitty, pearson"    , repel = TRUE)  )
  print(  fviz_dend(res.hc8_3three,   type="phylogenic", main="mcquitty, spearman"   , repel = TRUE)  )          
  dev.off() 
  
  # Compute hierarchical clustering by "median"
  res.hc1_3three <- hclust(res.dist1_3three, method = "median"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "median"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "median"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "median"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "median"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "median"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "median"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "median"  ) 
  
  pdf( file = paste(path_temp1, "8A_hierarchical-median-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="median, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="median, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="median, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="median, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="median, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="median, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="median, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="median, spearman"   )  )
  dev.off() 
  
  # Compute hierarchical clustering by "centroid"
  res.hc1_3three <- hclust(res.dist1_3three, method = "centroid"  )   
  res.hc2_3three <- hclust(res.dist2_3three, method = "centroid"  )   
  res.hc3_3three <- hclust(res.dist3_3three, method = "centroid"  )  
  res.hc4_3three <- hclust(res.dist4_3three, method = "centroid"  )   
  res.hc5_3three <- hclust(res.dist5_3three, method = "centroid"  )  
  res.hc6_3three <- hclust(res.dist6_3three, method = "centroid"  )  
  res.hc7_3three <- hclust(res.dist7_3three, method = "centroid"  )  
  res.hc8_3three <- hclust(res.dist8_3three, method = "centroid"  ) 
  
  pdf( file = paste(path_temp1, "9A_hierarchical-centroid-rectangle.pdf",  sep="/") )
  print(  fviz_dend(res.hc1_3three,   type="rectangle", main="centroid, euclidean"  )  )
  print(  fviz_dend(res.hc2_3three,   type="rectangle", main="centroid, maximum"    )  )
  print(  fviz_dend(res.hc3_3three,   type="rectangle", main="centroid, manhattan"  )  )
  print(  fviz_dend(res.hc4_3three,   type="rectangle", main="centroid, canberra"   )  )
  print(  fviz_dend(res.hc5_3three,   type="rectangle", main="centroid, binary"     )  )
  print(  fviz_dend(res.hc6_3three,   type="rectangle", main="centroid, minkowski"  )  )
  print(  fviz_dend(res.hc7_3three,   type="rectangle", main="centroid, pearson"    )  )
  print(  fviz_dend(res.hc8_3three,   type="rectangle", main="centroid, spearman"   )  )
  dev.off()  
}


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
MyPrcompObj_1_g <- function(  prcompObj2,   path2,   file2,  dataFrame_temp2   ) {
  
  sink( file = paste(path2, "/1A_PCA_results_",  file2,  ".txt",  sep="") )
  print( prcompObj2 )
  print( names(prcompObj2) )
  sink()
  
  sink( file = paste(path2, "/1B_PCA_summary_",  file2,  ".txt",   sep="") )
  print( summary(prcompObj2) )
  sink()
  
  sink( file = paste(path2, "/2_PCA_all_",  file2,  ".txt",   sep="") )
  print("####################### prcompObj2$sdev #########################")
  print(prcompObj2$sdev)
  print("####################### prcompObj2$rotation #########################")
  print(prcompObj2$rotation)
  print("####################### prcompObj2$center #########################")
  print(prcompObj2$center)
  print("####################### prcompObj2$scale #########################")
  print(prcompObj2$scale)
  print("####################### prcompObj2$x #########################")
  print(prcompObj2$x)
  sink()
  
  pdf( file = paste(path2, "/3_PCA_info_",  file2,  ".pdf",   sep="")  )
  plot(prcompObj2, type="lines")
  print( fviz_eig(prcompObj2) )
  dev.off() 
  
  my_fviz_pca_ind1 <- fviz_pca_ind(prcompObj2,
                                   col.ind = "cos2", # Color by the quality of representation
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   repel = TRUE     # Avoid text overlapping
  )
  my_fviz_pca_ind2 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE
  )
  my_fviz_pca_ind3 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none"
  )
  my_fviz_pca_ind4 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   ggtheme = MyTheme_1_g()
  )
  my_fviz_pca_ind5 <- fviz_pca_ind(prcompObj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none",
                                   ggtheme = MyTheme_1_g()
  )
  
  pdf( file=paste(path2, "/4_PCA-2D_", file2, ".pdf",  sep="") )
  print(my_fviz_pca_ind1)
  print(my_fviz_pca_ind2)
  print(my_fviz_pca_ind3)
  print(my_fviz_pca_ind4)
  print(my_fviz_pca_ind5)
  dev.off() 
  
  #############################
  prcompObj2_matrix <- prcompObj2$x
  prcompObj2_Contri  <- (prcompObj2$sdev)^2
  prcompObj2_Contri  <- prcompObj2_Contri/sum(prcompObj2_Contri)
  prcompObj2_Contri  <- prcompObj2_Contri * 100
  prcompObj2_Contri  <- round(prcompObj2_Contri, 2)
  
  label1_2two <-   paste( "PC1 ",  "(", prcompObj2_Contri[1], "%)", sep="" )
  label2_2two <-   paste( "PC2 ",  "(", prcompObj2_Contri[2], "%)", sep="" )
  label3_2two <-   paste( "PC3 ",  "(", prcompObj2_Contri[3], "%)", sep="" ) 
  
  dataframeA_2two  <- data.frame( as.data.frame(prcompObj2_matrix), mySex= as.vector(dataFrame_temp2$mysex), 
                                  myTech=as.vector(dataFrame_temp2$mytech),    myLabel=as.vector(dataFrame_temp2$mysampleID)   ) 
  
  FigureTemp1_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=3, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1_2two ,  path1=path2, fileName1="5A_PCA-PC1-PC2",  height1=2.5,  width1=4.3)
  
  FigureTemp2_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=5, alpha=0.5  )+ xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2_2two ,  path1=path2, fileName1="5B_PCA-PC1-PC2-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp4 <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=6, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +     
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4_2two ,  path1=path2, fileName1="5D_PCA-PC1-PC2-big",   height1=2.5,  width1=4.3)
  
  FigureTemp5_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=2, alpha=0.7  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5_2two ,  path1=path2, fileName1="5E_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  FigureTemp6_2two  <- ggplot( data = dataframeA_2two, aes(x = PC1, y = PC2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=5, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp6_2two ,  path1=path2, fileName1="5F_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  
  ## PC1, 2, and 3     dev.off() 
  pdf( file = paste(path2, "6_PCA-3d-by-scatter3D.pdf",  sep="/") )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),     
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),          
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ##############
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ######
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  dev.off()
}


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
MyPCAobj_FactoMineR_g <- function(  PCAobj2,   path2,   file2,  dataFrame_temp2   ) {
  
  sink( file = paste(path2, "/1A_PCA_results_",  file2,  ".txt",  sep="") )
  print( PCAobj2 )
  print( names(PCAobj2) )
  sink()
  
  sink( file = paste(path2, "/1B_PCA_summary_",  file2,  ".txt",   sep="") )
  print( summary(PCAobj2) )
  print( "#################################"  )
  myeig.val_PCAobj2 <- get_eigenvalue( PCAobj2 )
  print(myeig.val_PCAobj2)
  sink()
  
  sink( file = paste(path2, "/2_PCA_all_",  file2,  ".txt",   sep="") )
  print("####################### PCAobj2$ind$coord #########################")
  print(PCAobj2$ind$coord)
  print("####################### PCAobj2$ind #########################")
  print(PCAobj2$ind)
  sink()
  
  pdf( file = paste(path2, "/3_PCA_info_",  file2,  ".pdf",   sep="")  )
  print( fviz_eig(PCAobj2, addlabels = TRUE ) )
  print( fviz_screeplot(X=PCAobj2, choice = "variance", geom = "line",
                        barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                        ncp = 10, addlabels = TRUE) )
  print( fviz_screeplot(X=PCAobj2, choice = "eigenvalue", geom = "line",
                        barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                        ncp = 10, addlabels = TRUE) )
  print( fviz_screeplot(X=PCAobj2, choice = "variance", geom = "bar",
                        barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                        ncp = 10, addlabels = TRUE) )
  print( fviz_screeplot(X=PCAobj2, choice = "eigenvalue", geom = "bar",
                        barfill = "steelblue", barcolor = "steelblue", linecolor = "black",
                        ncp = 10, addlabels = TRUE) )
  dev.off() 
  
  
  my_fviz_pca_ind1 <- fviz_pca_ind(PCAobj2,
                                   col.ind = "cos2", # Color by the quality of representation
                                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                                   repel = TRUE     # Avoid text overlapping
  )
  my_fviz_pca_ind2 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE
  )
  my_fviz_pca_ind3 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none"
  )
  my_fviz_pca_ind4 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   ggtheme = MyTheme_1_g()
  )
  my_fviz_pca_ind5 <- fviz_pca_ind(PCAobj2,
                                   col.ind =  as.vector(as.factor( dataFrame_temp2$mytreatment) ), # color by groups
                                   repel = TRUE,
                                   label = "none",
                                   ggtheme = MyTheme_1_g()
  )
  
  pdf( file=paste(path2, "/4_PCA-2D_", file2, ".pdf",  sep="") )
  print(my_fviz_pca_ind1)
  print(my_fviz_pca_ind2)
  print(my_fviz_pca_ind3)
  print(my_fviz_pca_ind4)
  print(my_fviz_pca_ind5)
  dev.off() 
  
  
  #############################
  PCAobj2_matrix <- PCAobj2$ind$coord 
  PCAobj2_Contri  <- (PCAobj2$eig)[,2]
  PCAobj2_Contri  <- PCAobj2_Contri/sum(PCAobj2_Contri)
  PCAobj2_Contri  <- PCAobj2_Contri * 100
  PCAobj2_Contri  <- round(PCAobj2_Contri, 2)
  
  label1_2two <-   paste( "PC1 ",  "(", PCAobj2_Contri[1], "%)", sep="" )
  label2_2two <-   paste( "PC2 ",  "(", PCAobj2_Contri[2], "%)", sep="" )
  label3_2two <-   paste( "PC3 ",  "(", PCAobj2_Contri[3], "%)", sep="" ) 
  
  dataframeA_2two  <- data.frame( as.data.frame(PCAobj2_matrix), mySex= as.vector(dataFrame_temp2$mysex), 
                                  myTech=as.vector(dataFrame_temp2$mytech),    myLabel=as.vector(dataFrame_temp2$mysampleID)   )  
  
  FigureTemp1_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=3, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1_2two ,  path1=path2, fileName1="5A_PCA-PC1-PC2",  height1=2.5,  width1=4.3)
  
  FigureTemp2_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=5, alpha=0.5  )+ xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp2_2two ,  path1=path2, fileName1="5B_PCA-PC1-PC2-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=1  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot",   height1=2.5,  width1=4.3)
  
  FigureTemp3_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=2, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp3_2two ,  path1=path2, fileName1="5C_PCA-PC1-PC2-smallDot-alpha",   height1=2.5,  width1=4.3)
  
  FigureTemp4 <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech) )) + 
    geom_point(size=6, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +     
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) +  
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4_2two ,  path1=path2, fileName1="5D_PCA-PC1-PC2-big",   height1=2.5,  width1=4.3)
  
  FigureTemp5_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=2, alpha=0.7  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel)) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp5_2two ,  path1=path2, fileName1="5E_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  FigureTemp6_2two  <- ggplot( data = dataframeA_2two, aes(x = Dim.1, y = Dim.2, shape=as.factor(mySex), color=as.factor(myTech), label=myLabel )) + 
    geom_point(size=5, alpha=0.5  ) + xlab(label1_2two) +   ylab(label2_2two) +   
    scale_colour_manual(values=myTech_All_color_g) +  scale_shape_manual(values=mySex_All_shape_g) + 
    MyTheme_1_g( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) + geom_text_repel(aes(label = myLabel  )) +
    guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp6_2two ,  path1=path2, fileName1="5F_PCA-PC1-PC2-Labels",   height1=2.5,  width1=4.3)
  
  
  ## PC1, 2, and 3     dev.off() 
  pdf( file = paste(path2, "6_PCA-3d-by-scatter3D.pdf",  sep="/") )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),     
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),           
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),            
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),          
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ##############
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "g" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  
  ######
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=0, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=20, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=30, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=10,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=45, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=20,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=30,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  scatter3D(x = dataframeA_2two[,1],  y = dataframeA_2two[,2],  z = dataframeA_2two[,3], 
            col = MyTech_color_g( as.vector(dataFrame_temp2$mytech) )  ,  pch = MySex_Shape_g( as.vector(dataFrame_temp2$mysex) ),         bty = "b2" ,      
            colvar = NULL,  cex = 2, theta=60, phi=45,  xlab = label1_2two,  ylab = label2_2two,  zlab = label3_2two  )
  
  dev.off() 
}


MyPCA_1A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1    ) {   
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=TRUE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, "filterByQuantile_1", sep="/")
  path1_2 = paste(path1, "filterByQuantile_2", sep="/")
  path1_3 = paste(path1, "filterByQuantile_3", sep="/")
  path1_4 = paste(path1, "filterByQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
}


MyPCA_2A_g <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1,  dataFrame_temp1   ) {  
  if( ! file.exists(path1) ) { dir.create(path1, recursive = TRUE) }
  pdf( file=paste(path1, "/", file1, ".pdf", sep="") , width=width1, height=height1  )
  myObjTemp1 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp2 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=FALSE,  comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp3 <- PCASamples(mymeth1, screeplot=TRUE, scale=FALSE,  center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  myObjTemp4 <- PCASamples(mymeth1, screeplot=TRUE, scale=TRUE,   center=TRUE,   comp=c(1,2),  transpose=TRUE,  sd.filter=TRUE,   sd.threshold=sdThres1,  filterByQuantile=FALSE, obj.return=TRUE )
  dev.off()
  path1_1 = paste(path1, "filterNoQuantile_1", sep="/")
  path1_2 = paste(path1, "filterNoQuantile_2", sep="/")
  path1_3 = paste(path1, "filterNoQuantile_3", sep="/")
  path1_4 = paste(path1, "filterNoQuantile_4", sep="/")
  if( ! file.exists(path1_1) ) { dir.create(path1_1, recursive = TRUE) }
  if( ! file.exists(path1_2) ) { dir.create(path1_2, recursive = TRUE) }
  if( ! file.exists(path1_3) ) { dir.create(path1_3, recursive = TRUE) }
  if( ! file.exists(path1_4) ) { dir.create(path1_4, recursive = TRUE) }
  MyPrcompObj_1_g(  prcompObj2=myObjTemp1,   path2=path1_1,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp2,   path2=path1_2,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp3,   path2=path1_3,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
  MyPrcompObj_1_g(  prcompObj2=myObjTemp4,   path2=path1_4,   file2=file1,  dataFrame_temp2=dataFrame_temp1   )
}



MyPCA_3A_g <- function(  mymeth2 ,  path2,   file2, width2, height2,  dataFrame_temp2   )  {
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1A.kept100percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0   ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1B.kept90percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.1 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1C.kept80percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.2 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1D.kept70percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.3 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1E.kept60percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.4 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1F.kept50percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.5 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1G.kept40percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.6 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1H.kept30percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.7 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1I.kept20percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.8 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1J.kept10percent.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.9 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1K.kept5percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.95 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1L.kept4percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.96 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1M.kept3percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.97 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1N.kept2percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.98 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1O.kept1percent.pdf",   sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.99 , dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1P.kept0.5percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.995 ,dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_1A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".ByQuantile.1Q.kept0.1percent.pdf", sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.999 ,dataFrame_temp1=dataFrame_temp2   )                     
  
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2A.sd0.pdf",     sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0    ,  dataFrame_temp1=dataFrame_temp2   )                      
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2B.sd0.01.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.01 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2C.sd0.05.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.05 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2D.sd0.10.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.10 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2E.sd0.15.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.15 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2F.sd0.18.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.18 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2G.sd0.20.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.20 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2H.sd0.22.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.22 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2I.sd0.25.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.25 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2J.sd0.28.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.28 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2K.sd0.30.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.30 ,  dataFrame_temp1=dataFrame_temp2   )                     
  MyPCA_2A_g(  mymeth1 = mymeth2 ,   path1 = path2, file1 = paste(file2, ".NoQuantile.2L.sd0.35.pdf",  sep=""),   width1 = width2,   height1 = height2,   sdThres1 = 0.35 ,  dataFrame_temp1=dataFrame_temp2   )                       
}


##########################
## Annotating differentially methylated bases or regions
myRefSeqGenes_g = "/home/yp/AnnotationBED/hg38_RefSeq_Genes.bed"
myCpGIslands_g  = "/home/yp/AnnotationBED/hg38_CpG_islands.bed"
myRepeats_g     = "/home/yp/AnnotationBED/hg38_Repeats_rmsk.bed"
myImprintedRegions1_g = "/home/yp/AnnotationBED/67.Regions.PlosGenetics.ImprintedGenes.hg38.bed"
myImprintedRegions2_g = "/home/yp/AnnotationBED/75Regions.GR.hg38.bed"
myImprintedRegions3_g = "/home/yp/AnnotationBED/369Regions.GR.hg38.bed"
myImprintedRegions4_g = "/home/yp/AnnotationBED/merge1.imprintedRegions.hg38.bed"
myImprintedRegions5_g = "/home/yp/AnnotationBED/merge2.imprintedRegions.hg38.bed"

gene.obj_g     = readTranscriptFeatures(myRefSeqGenes_g)
cpg.obj_g      = readFeatureFlank(myCpGIslands_g, flank=10000, feature.flank.name=c("CpGi", "shores"))
myrepeat.obj_g = readFeatureFlank(myRepeats_g,    flank=10000, feature.flank.name=c("Repeats", "shores"))
imprint1.obj_g = readFeatureFlank(myImprintedRegions1_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint2.obj_g = readFeatureFlank(myImprintedRegions2_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint3.obj_g = readFeatureFlank(myImprintedRegions3_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint4.obj_g = readFeatureFlank(myImprintedRegions4_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
imprint5.obj_g = readFeatureFlank(myImprintedRegions5_g, flank=10000, feature.flank.name=c("ImprintedRegions", "shores"))
##########################


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myDiff_DMC_DMR_g <- function(  methobj2,   path2  ) {
  
  myDiff_2two_sub2 = calculateDiffMeth(methobj2,  num.cores=16)
  dim(myDiff_2two_sub2)
  names(myDiff_2two_sub2)
  head(myDiff_2two_sub2)
  
  myQvalue_2two_sub2 = myDiff_2two_sub2$qvalue
  myMethDi_2two_sub2 = myDiff_2two_sub2$meth.diff
  
  pdf(paste(path2,  "1A_qvalue_distribution.pdf",  sep="/"))
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1), freq=FALSE) )
  print( hist( myQvalue_2two_sub2[myQvalue_2two_sub2<0.1],  nclass=20, xlim=c(0, 0.1), freq=FALSE) )
  print( hist( myQvalue_2two_sub2,  nclass=100, xlim=c(0, 1), freq=TRUE) )
  print( hist( myQvalue_2two_sub2[myQvalue_2two_sub2<0.1],  nclass=20, xlim=c(0, 0.1), freq=TRUE) )
  dev.off() 
  
  pdf(paste(path2,  "1B_methDiff_distribution.pdf",  sep="/"))
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100), freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 50), freq=FALSE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=20], nclass=81, xlim=c(20, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=10], nclass=91, xlim=c(10, 100),  freq=FALSE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2<=10], nclass=100, xlim=c(0, 10),   freq=FALSE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 100), freq=TRUE) )
  print( hist(myMethDi_2two_sub2, nclass=100, xlim=c(0, 50), freq=TRUE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=20], nclass=81, xlim=c(20, 100),  freq=TRUE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2>=10], nclass=91, xlim=c(10, 100),  freq=TRUE) )
  print( hist(myMethDi_2two_sub2[myMethDi_2two_sub2<=10], nclass=100, xlim=c(0, 10),   freq=TRUE) )
  dev.off() 
  
  myColor1_2two_sub2 <- rep( "no",   times= length(myQvalue_2two_sub2) )
  myColor1_2two_sub2[ (abs(myMethDi_2two_sub2)>10) & (myQvalue_2two_sub2<0.001) ]  <- "yes"
  length( myColor1_2two_sub2[myColor1_2two_sub2=="yes"] )
  
  DataFrame2_2two_sub2 <- data.frame(myx1 = myMethDi_2two_sub2,   
                                     myy1 =  -log10(myQvalue_2two_sub2),  
                                     mycolor1 = myColor1_2two_sub2 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  
  ggsave( filename = paste(path2,  "1C_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  
  ggsave( filename = paste(path2,  "1C_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)  
  ggsave( filename = paste(path2,  "1C_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + ylim(0, 200)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-1.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-40, 40)  + ylim(0, 100)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-2.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  ggplot(DataFrame2_2two_sub2, aes(x=myx1, y=myy1, color= mycolor1 ) ) + 
    geom_point( shape = 20, alpha = 0.1 ) + scale_colour_manual(values=c("no"="black", "yes"="red")) +
    xlab( "Difference (%)" ) + ylab( "-log10(q-value)" ) + MyTheme_1_g(textSize1=14)  + xlim(-20, 20)   + ylim(0, 30)  
  ggsave( filename = paste(path2,  "1D_methDiff_qvalue-3.png",  sep="/"),  height=4, width=6, dpi = 1200 )
  
  myDiff25p.hypo_2two_sub2  = getMethylDiff(myDiff_2two_sub2, difference=10, qvalue=0.001, type="hypo" )  ## less enrich in ART
  myDiff25p.hyper_2two_sub2 = getMethylDiff(myDiff_2two_sub2, difference=10, qvalue=0.001, type="hyper")  ## more enrich in ART
  myDiff25p_2two_sub2       = getMethylDiff(myDiff_2two_sub2, difference=10, qvalue=0.05)
  myDiffTemp_2two_sub2      = getMethylDiff(myDiff_2two_sub2, difference=0,  qvalue=0.05)
  
  write.table(myDiff_2two_sub2 , file = paste(path2,"2A_diffMe-allsites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hypo_2two_sub2 , file = paste(path2,"2B_diffMe-hypo.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p.hyper_2two_sub2 , file = paste(path2,"2C_diffMe-hyper.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiff25p_2two_sub2 , file = paste(path2,"2D_AlldiffMesites.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  write.table(myDiffTemp_2two_sub2 , file = paste(path2,"2E_AlldiffMesites_q0.05_diff0.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  sink( file=paste(path2, "3A-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.txt", sep="/")   )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=FALSE,  qvalue.cutoff=0.001, meth.cutoff=10) )
  sink()
  
  pdf( file=paste(path2, "3B-distribution-of-hypoORhyper-methylated-bases-regions-per-chromosome.pdf", sep="/"), width=8, height=8    )
  print( diffMethPerChr(myDiff_2two_sub2,  plot=TRUE,  qvalue.cutoff=0.001, meth.cutoff=10) )
  dev.off()
  
  
  ########################## annotation for hypo sites.
  diffGeneAnn_2two_sub2_hypo = annotateWithGeneParts(as(myDiff25p.hypo_2two_sub2,  "GRanges"),  gene.obj_g)
  
  sink( file=paste(path2, "4A-distribution-onGenes-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_hypo,  percentage=TRUE)
  print(diffGeneAnn_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4B-distribution-onGenes-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffGeneAnn_2two_sub2_hypo,precedence=TRUE, main="differential methylation annotation") )
  dev.off()
  
  
  
  ##
  diffCpGann_2two_sub2_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2,"GRanges"),
                                                       cpg.obj_g$CpGi,  cpg.obj_g$shores,
                                                       feature.name="CpGi",flank.name="shores")
  
  sink( file=paste(path2, "4C-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffCpGann_2two_sub2_hypo,  percentage=TRUE)
  print(diffCpGann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4D-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffCpGann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffrepeatann_2two_sub2_hypo = annotateWithFeatureFlank(as(myDiff25p.hypo_2two_sub2,"GRanges"),
                                                          myrepeat.obj_g$Repeats,  myrepeat.obj_g$shores,
                                                          feature.name="Repeats",flank.name="shores")
  
  sink( file=paste(path2, "4E-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffrepeatann_2two_sub2_hypo,  percentage=TRUE)
  print(diffrepeatann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4F-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffrepeatann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted1Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5A-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5B-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted2Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5C-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5D-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted3Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5E-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5F-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted4Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5G-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5H-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted5Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=as(myDiff25p.hypo_2two_sub2, "GRanges"),
                                                              feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5I-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5J-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  ####################
  
  
  
  
  
  
  ########################## annotation for hyper sites.
  diffGeneAnn_2two_sub2_hyper = annotateWithGeneParts(as(myDiff25p.hyper_2two_sub2,"GRanges"),  gene.obj_g)
  
  sink( file=paste(path2, "10A-distribution-onGenes-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_hyper,  percentage=TRUE)
  print(diffGeneAnn_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "10B-distribution-onGenes-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffGeneAnn_2two_sub2_hyper,precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffCpGann_2two_sub2_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2,"GRanges"),
                                                        cpg.obj_g$CpGi,  cpg.obj_g$shores,
                                                        feature.name="CpGi",flank.name="shores")
  
  sink( file=paste(path2, "10C-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffCpGann_2two_sub2_hyper,  percentage=TRUE)
  print(diffCpGann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "10D-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffCpGann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffrepeatann_2two_sub2_hyper = annotateWithFeatureFlank(as(myDiff25p.hyper_2two_sub2,"GRanges"),
                                                           myrepeat.obj_g$Repeats,  myrepeat.obj_g$shores,
                                                           feature.name="Repeats",flank.name="shores")
  
  sink( file=paste(path2, "10E-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffrepeatann_2two_sub2_hyper,  percentage=TRUE)
  print(diffrepeatann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "10F-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffrepeatann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted1Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11A-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11B-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted2Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11C-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11D-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted3Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11E-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11F-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted4Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11G-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11H-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted5Ann_2two_sub2_hyper = annotateWithFeatureFlank(target=as(myDiff25p.hyper_2two_sub2, "GRanges"),
                                                               feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                               feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "11I-distribution-on-hyper.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hyper,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hyper)
  sink()
  
  pdf( file=paste(path2, "11J-distribution-onCpGs-hyper.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hyper, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
}


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction_1_g  <- function(  myobj_temp1,   path_temp1,   binSize_temp1, binBases_temp1, dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }
  myobj_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  tiles_2two = tileMethylCounts( myobj_2two,   win.size=binSize_temp1,   step.size=binSize_temp1,   cov.bases = binBases_temp1  )    
  meth_2two  = unite( tiles_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2_dimensions-tiles.txt", sep="/")  )
  print( tiles_2two )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf( file=paste(path_temp1_sub1, "4A_MethylationStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getMethylationStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4B_MethylationStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste(as.vector(dataFrame_temp1$mysampleID)[i],  ":", sep="") )
    print( getMethylationStats( tiles_2two[[i]] )  )
  }
  sink()
  
  pdf( file=paste(path_temp1_sub1, "5A_CoverageStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getCoverageStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  sink( file=paste(path_temp1_sub1, "5B_CoverageStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste( as.vector(dataFrame_temp1$mysampleID)[i],   ":", sep="") )
    print( getCoverageStats( tiles_2two[[i]] )  )
  }
  sink()
  
  #sink( file=paste(path_temp1_sub1, "6A_Correlation-tiles.txt", sep="/")  )
  #pdf( file=paste(path_temp1_sub1, "6A_Correlation-tiles.pdf", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=TRUE  )
  #getCorrelation(meth_2two, method = "spearman",  plot=TRUE  )
  #dev.off()
  #sink()
  
  sink( file=paste(path_temp1_sub1, "6B_pearsonCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "pearson",   plot=FALSE  )
  sink()
  
  sink( file=paste(path_temp1_sub1, "6C_spearmanCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "spearman",  plot=FALSE  )
  sink()
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=8,   height2=5 )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=4,   height2=4 )
  
  #path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  #if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  #MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=4,   height2=4,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   dataFrame_temp1=dataFrame_temp1  )
  
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )
  
  
} 


#dataFrame_temp111 <- data.frame(
#  mysampleID  = c(mySampleID_NC_g,  mySampleID_IVF_fresh_g),
#  mytreatment = c(myTreatment_NC_g, myTreatment_IVF_fresh_g),
#  mysex       = c(Sex_NC_g,  Sex_IVF_fresh_g),
#  mytech      = c(Tech_NC_g, Tech_IVF_fresh_g)    
#)
myMainFunction_2_g  <- function(  myobj_temp1,   path_temp1,   binSize_temp1, binBases_temp1, dataFrame_temp1  ) {
  if( ! file.exists(path_temp1) ) { dir.create(path_temp1, recursive = TRUE) }
  myobj_2two <- reorganize(myobj_temp1,  sample.ids = as.vector(dataFrame_temp1$mysampleID),  treatment  = as.vector(dataFrame_temp1$mytreatment) )
  
  path_temp1_sub1 = paste(path_temp1, "1_stats_information", sep="/")
  if( ! file.exists(path_temp1_sub1) ) { dir.create(path_temp1_sub1, recursive = TRUE) }
  
  sink( file=paste(path_temp1_sub1,  "1_select-subSets.txt", sep="/")  )
  print( length(myobj_2two) )
  print( getSampleID(myobj_2two) )
  print( getTreatment(myobj_2two) )
  sink()
  
  tiles_2two = tileMethylCounts( myobj_2two,   win.size=binSize_temp1,   step.size=binSize_temp1,   cov.bases = binBases_temp1  )    
  meth_2two  = unite( tiles_2two, destrand=FALSE, mc.cores=16   )   ## 100% overlap
  mat_2two   = percMethylation( meth_2two )
  
  sink( file=paste(path_temp1_sub1 , "2_dimensions-tiles.txt", sep="/")  )
  print( tiles_2two )
  print("#########dimensions:")
  print( dim(meth_2two)  )   
  print( dim(mat_2two)   )
  sink()
  
  write.table(meth_2two , 
              file = paste(path_temp1_sub1,   "3A_meth-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  write.table(mat_2two , 
              file = paste(path_temp1_sub1,   "3B_mat-tiles.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf( file=paste(path_temp1_sub1, "4A_MethylationStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getMethylationStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  
  sink( file=paste(path_temp1_sub1, "4B_MethylationStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste(as.vector(dataFrame_temp1$mysampleID)[i],  ":", sep="") )
    print( getMethylationStats( tiles_2two[[i]] )  )
  }
  sink()
  
  pdf( file=paste(path_temp1_sub1, "5A_CoverageStats-tiles.pdf", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    getCoverageStats(tiles_2two[[i]], plot=TRUE, both.strands=FALSE )
  }
  dev.off()
  sink( file=paste(path_temp1_sub1, "5B_CoverageStats-tiles.txt", sep="/")  )
  for( i in c(1:length(myobj_2two)) ) {
    print("##############################################")
    print( paste( as.vector(dataFrame_temp1$mysampleID)[i],   ":", sep="") )
    print( getCoverageStats( tiles_2two[[i]] )  )
  }
  sink()
  
  #sink( file=paste(path_temp1_sub1, "6A_Correlation-tiles.txt", sep="/")  )
  #pdf( file=paste(path_temp1_sub1, "6A_Correlation-tiles.pdf", sep="/")  )
  #getCorrelation(meth_2two, method = "pearson",   plot=TRUE  )
  #getCorrelation(meth_2two, method = "spearman",  plot=TRUE  )
  #dev.off()
  #sink()
  
  sink( file=paste(path_temp1_sub1, "6B_pearsonCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "pearson",   plot=FALSE  )
  sink()
  
  sink( file=paste(path_temp1_sub1, "6C_spearmanCorrelation-tiles.txt", sep="/")  )
  getCorrelation(meth_2two, method = "spearman",  plot=FALSE  )
  sink()
  
  
  path_temp1_sub2 = paste(path_temp1, "2_HierarchicalClustering_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub2) ) { dir.create(path_temp1_sub2, recursive = TRUE) }
  MyCluster_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub2,     file2="HierarchicalClustering_byMethylKit_",   width2=8,   height2=5 )
  
  path_temp1_sub3 = paste(path_temp1, "3_PCA_byMethylKit", sep="/")
  if( ! file.exists(path_temp1_sub3) ) { dir.create(path_temp1_sub3, recursive = TRUE) }
  MyPCA_3_g(  mymeth2=meth_2two ,  path2=path_temp1_sub3,     file2="PCA_byMethylKit_",   width2=4,   height2=4 )
  
  #path_temp1_sub4 = paste(path_temp1, "4_PCAinfor_byMethylKit", sep="/")
  #if( ! file.exists(path_temp1_sub4) ) { dir.create(path_temp1_sub4, recursive = TRUE) }
  #MyPCA_3A_g(  mymeth2=meth_2two ,  path2=path_temp1_sub4,     file2="PCAinfor_byMethylKit_",   width2=4,   height2=4,  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub5 = paste(path_temp1, "5_PCA_DirectlyByPrcomp", sep="/")
  if( ! file.exists(path_temp1_sub5) ) { dir.create(path_temp1_sub5, recursive = TRUE) }
  PCA_2two_sub5 <- prcomp( t(mat_2two) )
  MyPrcompObj_1_g(  prcompObj2=PCA_2two_sub5,   path2=path_temp1_sub5,   file2="PCA_DirectlyByPrcomp",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub6 = paste(path_temp1, "6_PCA_byFactoMineR", sep="/")
  if( ! file.exists(path_temp1_sub6) ) { dir.create(path_temp1_sub6, recursive = TRUE) }
  PCA_2two_sub6 <- PCA( t(mat_2two) , graph=FALSE)
  MyPCAobj_FactoMineR_g(  PCAobj2=PCA_2two_sub6,   path2=path_temp1_sub6,   file2="PCA_byFactoMineR",  dataFrame_temp2=dataFrame_temp1   )
  
  path_temp1_sub7 = paste(path_temp1, "7_HierarchicalClustering", sep="/")
  if( ! file.exists(path_temp1_sub7) ) { dir.create(path_temp1_sub7, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=mat_2two,   path_temp1=path_temp1_sub7,   dataFrame_temp1=dataFrame_temp1  )
  
  path_temp1_sub8 = paste(path_temp1, "8_DMR", sep="/")
  if( ! file.exists(path_temp1_sub8) ) { dir.create(path_temp1_sub8, recursive = TRUE) }
  myDiff_DMC_DMR_g(  methobj2=meth_2two,   path2=path_temp1_sub8  )
  
  
} 



## scatter Diagram for 2 samples and colour each category.
MyScatterDiagram_1 <- function(vector2X,  vector2Y,  path2,   fileName2,   xLab2,   yLab2,  title2,  height2=4,  width2=4,  yMin2=0, yMax2=2,  xMin2=0, xMax2=2,  alpha2=0.5, diffThres2=1.5, colours2=c("red", "blue", "purple")) {                     
  vector2Y[vector2Y>yMax2] <- yMax2
  vector2Y[vector2Y<yMin2] <- yMin2
  vector2X[vector2X>xMax2] <- xMax2
  vector2X[vector2X<xMin2] <- xMin2
  myRatio  <- (vector2X - vector2Y) 
  myRatio1 <- myRatio
  myRatio2 <- myRatio
  myRatio3 <- myRatio
  myRatio4 <- myRatio
  myRatio5 <- myRatio
  
  myRatio1[myRatio1 >= diffThres2] <- 200 
  myRatio1[myRatio1 <= -diffThres2] <- -200 
  myRatio1[(myRatio1 < diffThres2)&(myRatio1 > -diffThres2)] <- 0
  print(myRatio1)
  
  dataframeA <- data.frame( xAxis = vector2X,  yAxis = vector2Y , SampleType = factor(myRatio1) ) 
  FigureTemp1 <- ggplot( data = dataframeA, aes(x = xAxis, y = yAxis,  color=factor(SampleType) )) + 
    geom_point(size=0.1, alpha=alpha2, shape=20 ) +  geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + scale_colour_manual( values=colours2 ) +
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_1_diffThres2",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  myRatio2[myRatio2 >= 50] <- 200
  myRatio2[myRatio2 <= -50] <- -200
  myRatio2[(myRatio2 < 50)&(myRatio2 > -50)] <- 0
  dataframeB <- data.frame( xAxis = vector2X,  yAxis = vector2Y , SampleType = factor(myRatio2) ) 
  FigureTemp1 <- ggplot( data = dataframeB, aes(x = xAxis, y = yAxis,  color=factor(SampleType) )) + 
    geom_point(size=0.1, alpha=alpha2, shape=20 ) +  geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + scale_colour_manual( values=colours2 ) +
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_2_diff50",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  myRatio3[myRatio3 >= 25] <- 200
  myRatio3[myRatio3 <= -25] <- -200
  myRatio3[(myRatio3 < 25)&(myRatio3 > -25)] <- 0
  dataframeC <- data.frame( xAxis = vector2X,  yAxis = vector2Y , SampleType = factor(myRatio3) ) 
  FigureTemp1 <- ggplot( data = dataframeC, aes(x = xAxis, y = yAxis,  color=factor(SampleType) )) + 
    geom_point(size=0.1, alpha=alpha2, shape=20 ) +  geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + scale_colour_manual( values=colours2 ) +
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_3_diff25",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  
  myRatio4[myRatio4 >= 20] <- 200
  myRatio4[myRatio4 <= -20] <- -200
  myRatio4[(myRatio4 < 20)&(myRatio4 > -20)] <- 0
  dataframeC <- data.frame( xAxis = vector2X,  yAxis = vector2Y , SampleType = factor(myRatio4) ) 
  FigureTemp1 <- ggplot( data = dataframeC, aes(x = xAxis, y = yAxis,  color=factor(SampleType) )) + 
    geom_point(size=0.1, alpha=alpha2, shape=20 ) +  geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + scale_colour_manual( values=colours2 ) +
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_4_diff20",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  
  myRatio5[myRatio5 >= 5] <- 200
  myRatio5[myRatio5 <= -5] <- -200
  myRatio5[(myRatio5 < 5)&(myRatio5 > -5)] <- 0
  dataframeC <- data.frame( xAxis = vector2X,  yAxis = vector2Y , SampleType = factor(myRatio5) ) 
  FigureTemp1 <- ggplot( data = dataframeC, aes(x = xAxis, y = yAxis,  color=factor(SampleType) )) + 
    geom_point(size=0.1, alpha=alpha2, shape=20 ) +  geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + scale_colour_manual( values=colours2 ) +
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_5_diff5",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
}



## 2D density plot
MyScatterDiagram_2 <- function(vector2X,  vector2Y,  path2,   fileName2,   xLab2,   yLab2,  title2,  height2=4,  width2=4,  yMin2=0, yMax2=2,  xMin2=0, xMax2=2 ) {
  vector2Y[vector2Y>yMax2] <- yMax2
  vector2Y[vector2Y<yMin2] <- yMin2
  vector2X[vector2X>xMax2] <- xMax2
  vector2X[vector2X<xMin2] <- xMin2
  
  dataframeA <- data.frame( xAxis = vector2X,  yAxis = vector2Y  )
  
  FigureTemp1 <- ggplot( data = dataframeA   ) +  
    stat_density_2d(  aes(x=xAxis, y=yAxis, fill = ..density.. ) ,  geom = "raster" , contour=FALSE ) + 
    geom_abline(slope=1, intercept=10, lty=2, col="black", size=0.3) + geom_abline(slope=1, intercept= -10, lty=2, col="black", size=0.3) +
    geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +  scale_fill_continuous( low=c("white", "blue") ,  high=c("yellow", "red") , na.value = "white" ) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_1_4colours",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp2 <- ggplot( data = dataframeA   ) +  
    stat_density_2d(  aes(x=xAxis, y=yAxis, fill = ..density.. ) ,  geom = "raster" , contour=FALSE ) + 
    geom_abline(slope=1, intercept=10, lty=2, col="black", size=0.3) + geom_abline(slope=1, intercept= -10, lty=2, col="black", size=0.3) +
    geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +  scale_fill_continuous( low=c("white", "blue") ,  high="red" , na.value = "white" ) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_2_3colours",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp3 <- ggplot( data = dataframeA   ) +  
    stat_density_2d(  aes(x=xAxis, y=yAxis, fill = ..density.. ) ,  geom = "raster" , contour=FALSE  ) + 
    geom_abline(slope=1, intercept=10, lty=2, col="black", size=0.3) + geom_abline(slope=1, intercept= -10, lty=2, col="black", size=0.3) +
    geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +  scale_fill_continuous( low= "white"  ,  high="red" , na.value = "white" ) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "_3_2colours",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp4 <- ggplot( data = dataframeA   ) +  
    stat_density_2d(  aes(x=xAxis, y=yAxis, fill = ..density.. ) ,  geom = "raster" , contour=FALSE  ) + 
    geom_abline(slope=1, intercept=10, lty=2, col="black", size=0.3) + geom_abline(slope=1, intercept= -10, lty=2, col="black", size=0.3) +
    geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +  scale_fill_continuous( low= c("white", "red")  ,  high="red4" , na.value = "white" ) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "_4_other3colours",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp5 <- ggplot( data = dataframeA   ) +  
    stat_density_2d(  aes(x=xAxis, y=yAxis, fill = ..density.. ) ,  geom = "raster" , contour=FALSE  ) + 
    geom_abline(slope=1, intercept=10, lty=2, col="black", size=0.3) + geom_abline(slope=1, intercept= -10, lty=2, col="black", size=0.3) +
    geom_abline(slope=1, intercept=0, lty=2, col="black", size=0.3) +  scale_fill_continuous( low= c("white", "blue")  ,  high=c("red", "red4") , na.value = "white" ) +
    ylim(yMin2, yMax2) +  xlim(xMin2, xMax2)+  xlab(xLab2) +   ylab(yLab2) +   ggtitle(title2) + 
    guides( colour = guide_legend(override.aes = list(size=10)) )  + MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL)
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_5_other4colours",   sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
} 


MyBoxViolinPlot_1 <- function(vector2,   sampleType2,  colours2,   path2,   fileName2,  title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) +  
    geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
    geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=NA,  outlier.size=0, size=0.5, fill=colours2 ) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray" ) +   
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 2) +  
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-2Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-3Adjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 4) +  
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.shape=NA,  outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="yellow4", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-4Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)  
}  


MyBoxViolinPlot_2 <- function(vector2,   sampleType2,  sampleType3,  colours2,   path2,   fileName2,  title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2,  sampleTypeB=sampleType3   ) 
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleTypeB) ) +  
    #geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
    geom_boxplot( outlier.shape=NA, outlier.size=0, notch=TRUE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis , fill=sampleTypeB) ) +  
    geom_violin( colour = NA  ) + 
    #geom_boxplot( outlier.shape=NA, outlier.size=0, size=1,   alpha=0.001   ) +   
    #stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_violinPlot",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp3 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis , fill=sampleTypeB) ) +  
    geom_violin(  colour = NA  ) + 
    #geom_boxplot( outlier.shape=NA, outlier.size=0, size=1,   alpha=0.001   ) +   
    stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "_violinPlot-mean",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp4 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis , fill=sampleTypeB) ) +  
    geom_violin(  colour = NA  ) + 
    geom_boxplot( outlier.shape=NA, outlier.size=0, size=1,   alpha=0.001   ) +   
    stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "_violinBoxPlot",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp5 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis , fill=sampleTypeB) ) +  
    geom_violin(colour = NA ,  adjust = 2 ) + 
    geom_boxplot( outlier.shape=NA, outlier.size=0, size=1,   alpha=0.001   ) +   
    stat_summary( position=position_dodge(width=0.75), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_violinBoxPlot-2Ajust",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
}  


myAnnotation_1 <- function( includeRegions2, path2 ) {  
  if( ! file.exists(path2) ) { dir.create(path2, recursive = TRUE) }
  
  Regions2_temp1 <- as(includeRegions2,  "GRanges")
  
  diffGeneAnn_2two_sub2_hypo = annotateWithGeneParts(Regions2_temp1,  gene.obj_g)
  
  sink( file=paste(path2, "1A-distribution-onGenes.txt", sep="/")   )
  getFeatsWithTargetsStats(diffGeneAnn_2two_sub2_hypo,  percentage=TRUE)
  print(diffGeneAnn_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "1B-distribution-onGenes.pdf", sep="/")   )
  print(plotTargetAnnotation(diffGeneAnn_2two_sub2_hypo,precedence=TRUE, main="differential methylation annotation") )
  dev.off()
  
  
  
  ##
  diffCpGann_2two_sub2_hypo = annotateWithFeatureFlank(Regions2_temp1,
                                                       cpg.obj_g$CpGi,  cpg.obj_g$shores,
                                                       feature.name="CpGi",flank.name="shores")
  
  sink( file=paste(path2, "2A-distribution-on-hypo.txt", sep="/")   )
  getFeatsWithTargetsStats(diffCpGann_2two_sub2_hypo,  percentage=TRUE)
  print(diffCpGann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "2B-distribution-onCpGs-hypo.pdf", sep="/")   )
  print(plotTargetAnnotation(diffCpGann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffrepeatann_2two_sub2_hypo = annotateWithFeatureFlank(Regions2_temp1,
                                                          myrepeat.obj_g$Repeats,  myrepeat.obj_g$shores,
                                                          feature.name="Repeats",flank.name="shores")
  
  sink( file=paste(path2, "3A-distribution-on-Repeats.txt", sep="/")   )
  getFeatsWithTargetsStats(diffrepeatann_2two_sub2_hypo,  percentage=TRUE)
  print(diffrepeatann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "3B-distribution-onRepeats.pdf", sep="/")   )
  print(plotTargetAnnotation(diffrepeatann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted1Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=Regions2_temp1,
                                                              feature=imprint1.obj_g$ImprintedRegions,  flank=imprint1.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "4A-distribution-on-ImprintedRegions.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted1Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted1Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "4B-distribution-onImprintedRegions.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted1Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted2Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=Regions2_temp1,
                                                              feature=imprint2.obj_g$ImprintedRegions,  flank=imprint2.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "5A-distribution-on-ImprintedRegions.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted2Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted2Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "5B-distribution-onImprintedRegions.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted2Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted3Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=Regions2_temp1,
                                                              feature=imprint3.obj_g$ImprintedRegions,  flank=imprint3.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "6A-distribution-on-ImprintedRegions.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted3Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted3Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "6B-distribution-onImprintedRegions.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted3Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  ##
  diffImprinted4Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=Regions2_temp1,
                                                              feature=imprint4.obj_g$ImprintedRegions,  flank=imprint4.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "7A-distribution-on-ImprintedRegions.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted4Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted4Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "7B-distribution-onImprintedRegions.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted4Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  
  
  
  
  ##
  diffImprinted5Ann_2two_sub2_hypo = annotateWithFeatureFlank(target=Regions2_temp1,
                                                              feature=imprint5.obj_g$ImprintedRegions,  flank=imprint5.obj_g$shores,
                                                              feature.name="ImprintedRegions",  flank.name="shores")
  
  sink( file=paste(path2, "8A-distribution-on-ImprintedRegions.txt", sep="/")   )
  getFeatsWithTargetsStats(diffImprinted5Ann_2two_sub2_hypo,  percentage=TRUE)
  print(diffImprinted5Ann_2two_sub2_hypo)
  sink()
  
  pdf( file=paste(path2, "8B-distribution-onImprintedRegions.pdf", sep="/")   )
  print(plotTargetAnnotation(diffImprinted5Ann_2two_sub2_hypo, precedence=TRUE, main="differential methylation annotation"))
  dev.off()
  ####################
  
}  



##################################################################################################################







## Define the sample groups about children
##################################################################################################################
Files_NC_children_g <- c(
  paste( inputDir_g,   "67_E24C-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "67_E24D-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56C-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56D-boy-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "69_NC-E123-C-Boy_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "69_NC-E123-D-Boy_Rep1.bismark.cov",  sep="/" ) 
)

mySampleID_NC_children_g <- rep( x="NC", times=length(Files_NC_children_g) )
Tech_NC_children_g <- mySampleID_NC_children_g

for(i in c(1:length(mySampleID_NC_children_g)) ) {
  mySampleID_NC_children_g[i] = paste(mySampleID_NC_children_g[i], i, sep="_")
}

myTreatment_NC_children_g <- rep( x=0,  times=length(Files_NC_children_g) )

Sex_NC_children_g = rep( x="boy",  times=length(Files_NC_children_g) )  
for(i in c(1:length(Sex_NC_children_g)) ) {
  Sex_NC_children_g[i] = "boy"
  if(i>=4) { Sex_NC_children_g[i] = "boy" }
}



Files_IVF_fresh_children_g <- c(
  paste( inputDir_g,   "70_E113C-boy-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "70_E113D-boy-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-C-Boy_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-D-Boy_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "72_ART-W779-C-Boy_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "72_W779D-ART-boy_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_IVF_fresh_children_g <- rep( x="IVF_fresh", times=length(Files_IVF_fresh_children_g) )
Tech_IVF_fresh_children_g <- mySampleID_IVF_fresh_children_g

for(i in c(1:length(mySampleID_IVF_fresh_children_g)) ) {
  mySampleID_IVF_fresh_children_g[i] = paste(mySampleID_IVF_fresh_children_g[i], i, sep="_")
}

myTreatment_IVF_fresh_children_g <- rep( x=1,  times=length(Files_IVF_fresh_children_g) )

Sex_IVF_fresh_children_g = rep( x="boy",      times=length(Files_IVF_fresh_children_g) ) 
for(i in c(1:length(Sex_IVF_fresh_children_g)) ) {
  Sex_IVF_fresh_children_g[i] = "boy"
  if(i>=4) { Sex_IVF_fresh_children_g[i] = "boy" }
}



Files_IVF_frozen_children_g <- c(
  paste( inputDir_g,   "9_W1365C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "9_Q5-W1365D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733C-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "11_W1398C-boy-IVF-frozen_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "11_W1398D-boy-IVF-frozen_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_IVF_frozen_children_g <- rep( x="IVF_frozen", times=length(Files_IVF_frozen_children_g) )
Tech_IVF_frozen_children_g <- mySampleID_IVF_frozen_children_g

for(i in c(1:length(mySampleID_IVF_frozen_children_g)) ) {
  mySampleID_IVF_frozen_children_g[i] = paste(mySampleID_IVF_frozen_children_g[i], i, sep="_")
}

myTreatment_IVF_frozen_children_g <- rep( x=2,  times=length(Files_IVF_frozen_children_g) )

Sex_IVF_frozen_children_g = rep( x="boy",      times=length(Files_IVF_frozen_children_g) )
for(i in c(1:length(Sex_IVF_frozen_children_g)) ) {
  Sex_IVF_frozen_children_g[i] = "boy"
  if(i>=4) { Sex_IVF_frozen_children_g[i] = "boy" }
}



Files_ICSI_fresh_children_g <- c(
  paste( inputDir_g,   "12_W1579C-boy-ICSI-fresh-merge_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "12_Q17-W1579D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647C-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "14_W1719C-boy-ICSI-fresh_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "14_W1719D-boy-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_fresh_children_g <- rep( x="ICSI_fresh", times=length(Files_ICSI_fresh_children_g) )
Tech_ICSI_fresh_children_g <- mySampleID_ICSI_fresh_children_g

for(i in c(1:length(mySampleID_ICSI_fresh_children_g)) ) {
  mySampleID_ICSI_fresh_children_g[i] = paste(mySampleID_ICSI_fresh_children_g[i], i, sep="_")
}

myTreatment_ICSI_fresh_children_g <- rep( x=3,  times=length(Files_ICSI_fresh_children_g) )

Sex_ICSI_fresh_children_g = rep( x="boy",      times=length(Files_ICSI_fresh_children_g) )  
for(i in c(1:length(Sex_ICSI_fresh_children_g)) ) {
  Sex_ICSI_fresh_children_g[i] = "boy"
  if(i>=4) { Sex_ICSI_fresh_children_g[i] = "boy" }
}



Files_ICSI_frozen_children_g <- c(
  paste( inputDir_g,   "15_Q6-W871C-boy-ICSI-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "15_Q4-W871D-boy-ICSI-frozen_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_frozen_children_g <- rep( x="ICSI_frozen", times=length(Files_ICSI_frozen_children_g) )
Tech_ICSI_frozen_children_g <- mySampleID_ICSI_frozen_children_g

for(i in c(1:length(mySampleID_ICSI_frozen_children_g)) ) {
  mySampleID_ICSI_frozen_children_g[i] = paste(mySampleID_ICSI_frozen_children_g[i], i, sep="_")
}

myTreatment_ICSI_frozen_children_g <- rep( x=4,  times=length(Files_ICSI_frozen_children_g) )

Sex_ICSI_frozen_children_g = rep( x="boy",      times=length(Files_ICSI_frozen_children_g) )  
for(i in c(1:length(Sex_ICSI_frozen_children_g)) ) {
  Sex_ICSI_frozen_children_g[i] = "boy"
  if(i>=4) { Sex_ICSI_frozen_children_g[i] = "boy" }
}





Files_All_vector_children_g <- c(
  Files_NC_children_g,
  Files_IVF_fresh_children_g,
  Files_IVF_frozen_children_g,
  Files_ICSI_fresh_children_g,
  Files_ICSI_frozen_children_g  
)
Files_All_list_children_g <- as.list( Files_All_vector_children_g )


mySampleID_All_vector_children_g <- c( 
  mySampleID_NC_children_g,
  mySampleID_IVF_fresh_children_g,
  mySampleID_IVF_frozen_children_g,
  mySampleID_ICSI_fresh_children_g,
  mySampleID_ICSI_frozen_children_g  
)
mySampleID_All_list_children_g <- as.list( mySampleID_All_vector_children_g )


myTreatment_All_vector_children_g <- c( 
  myTreatment_NC_children_g,
  myTreatment_IVF_fresh_children_g,
  myTreatment_IVF_frozen_children_g,
  myTreatment_ICSI_fresh_children_g,
  myTreatment_ICSI_frozen_children_g 
)       
myTreatment_All_list_children_g <- as.list( myTreatment_All_vector_children_g )


mySex_All_vector_children_g <- c( 
  Sex_NC_children_g,
  Sex_IVF_fresh_children_g,
  Sex_IVF_frozen_children_g,
  Sex_ICSI_fresh_children_g,
  Sex_ICSI_frozen_children_g 
)       
mySex_All_list_children_g <- as.list( mySex_All_vector_children_g )


myTech_All_vector_children_g <- c( 
  Tech_NC_children_g,
  Tech_IVF_fresh_children_g,
  Tech_IVF_frozen_children_g,
  Tech_ICSI_fresh_children_g,
  Tech_ICSI_frozen_children_g 
)       
myTech_All_list_children_g <- as.list( myTech_All_vector_children_g )


## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple
mySex_All_shape_children_g  = c( "boy"=16, "girl"=17, "father"=1, "mother"=2 ) 
myTech_All_color_children_g = c( "NC"="black", "IVF_fresh"="blue", "IVF_frozen"="green",  "ICSI_fresh"="red", "ICSI_frozen"="purple" )
mySex_All_shape_children_g
myTech_All_color_children_g

MySex_Shape_children_g <- function(  mySex_vector   ) {
  mySex_shape2  = mySex_vector
  for(i in c(1:length(mySex_shape2)) ) {
    if(mySex_shape2[i] == "boy")    { mySex_shape2[i] = c("boy"=16)   }
    if(mySex_shape2[i] == "girl")   { mySex_shape2[i] = c("girl"=17)  }
    if(mySex_shape2[i] == "father") { mySex_shape2[i] = c("father"=1) }
    if(mySex_shape2[i] == "mother") { mySex_shape2[i] = c("mother"=2) }
  }
  mySex_shape2 = as.numeric(mySex_shape2)
  names(mySex_shape2) = mySex_vector
  return(mySex_shape2)
}

MyTech_color_children_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == "NC")          { myTech_color[i] = c("NC"="black") }
    if(myTech_color[i] == "IVF_fresh")   { myTech_color[i] = c("IVF_fresh"="blue") }
    if(myTech_color[i] == "IVF_frozen")  { myTech_color[i] = c("IVF_frozen"="green") }
    if(myTech_color[i] == "ICSI_fresh")  { myTech_color[i] = c("ICSI_fresh"="red") }
    if(myTech_color[i] == "ICSI_frozen") { myTech_color[i] = c("ICSI_frozen"="purple") }
  }
  names(myTech_color) = myTech_vector
  return(myTech_color)
}

MySex_Shape_children_g( mySex_All_vector_children_g )
MyTech_color_children_g( myTech_All_vector_children_g )
## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple


length( Files_All_vector_children_g )
length( Files_All_list_children_g )
length( mySampleID_All_vector_children_g )
length( mySampleID_All_list_children_g )
length( myTreatment_All_vector_children_g )
length( myTreatment_All_list_children_g )
length( mySex_All_vector_children_g )
length( mySex_All_list_children_g )
length( myTech_All_vector_children_g )
length( myTech_All_list_children_g )
##################################################################################################################







## Define the sample groups about parents
##################################################################################################################

Files_NC_parents_g <- c(
  paste( inputDir_g,   "67_NC-E24-F-Father_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56F-father-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "69_NC-E123-F-Father_Rep1.bismark.cov",  sep="/" ),
  paste( inputDir_g,   "67_NC-E24-M-Mother_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "68_E56M-mother-NC_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "69_NC-E123-M-Mother-merge_Rep1.bismark.cov",  sep="/" ) 
)

mySampleID_NC_parents_g <- rep( x="NC", times=length(Files_NC_parents_g) )
Tech_NC_parents_g <- mySampleID_NC_parents_g

for(i in c(1:length(mySampleID_NC_parents_g)) ) {
  mySampleID_NC_parents_g[i] = paste(mySampleID_NC_parents_g[i], i, sep="_")
}

myTreatment_NC_parents_g <- rep( x=0,  times=length(Files_NC_parents_g) )

Sex_NC_parents_g = rep( x="father",  times=length(Files_NC_parents_g) )  
for(i in c(1:length(Sex_NC_parents_g)) ) {
  Sex_NC_parents_g[i] = "father"
  if(i>=4) { Sex_NC_parents_g[i] = "mother" }
}



Files_IVF_fresh_parents_g <- c(
  paste( inputDir_g,   "70_E113F-father-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-F-Father_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "72_ART-W779-F-Father_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "70_E113M-mother-ART_Rep3.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "71_ART-W58-M-Mother_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "72_ART-W779-M-Mother_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_IVF_fresh_parents_g <- rep( x="IVF_fresh", times=length(Files_IVF_fresh_parents_g) )
Tech_IVF_fresh_parents_g <- mySampleID_IVF_fresh_parents_g

for(i in c(1:length(mySampleID_IVF_fresh_parents_g)) ) {
  mySampleID_IVF_fresh_parents_g[i] = paste(mySampleID_IVF_fresh_parents_g[i], i, sep="_")
}

myTreatment_IVF_fresh_parents_g <- rep( x=1,  times=length(Files_IVF_fresh_parents_g) )

Sex_IVF_fresh_parents_g = rep( x="father",      times=length(Files_IVF_fresh_parents_g) ) 
for(i in c(1:length(Sex_IVF_fresh_parents_g)) ) {
  Sex_IVF_fresh_parents_g[i] = "father"
  if(i>=4) { Sex_IVF_fresh_parents_g[i] = "mother" }
}



Files_IVF_frozen_parents_g <- c(
  paste( inputDir_g,   "9_W1365F-father-IVF-frozen-merge_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733F-father-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "11_W1398F-father-IVF-frozen_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "9_Q1-W1365M-mother-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "10_W1733M-Mother-IVF-frozen_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "12_Q19-W1579M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_IVF_frozen_parents_g <- rep( x="IVF_frozen", times=length(Files_IVF_frozen_parents_g) )
Tech_IVF_frozen_parents_g <- mySampleID_IVF_frozen_parents_g

for(i in c(1:length(mySampleID_IVF_frozen_parents_g)) ) {
  mySampleID_IVF_frozen_parents_g[i] = paste(mySampleID_IVF_frozen_parents_g[i], i, sep="_")
}

myTreatment_IVF_frozen_parents_g <- rep( x=2,  times=length(Files_IVF_frozen_parents_g) )

Sex_IVF_frozen_parents_g = rep( x="father",      times=length(Files_IVF_frozen_parents_g) )
for(i in c(1:length(Sex_IVF_frozen_parents_g)) ) {
  Sex_IVF_frozen_parents_g[i] = "father"
  if(i>=4) { Sex_IVF_frozen_parents_g[i] = "mother" }
}



Files_ICSI_fresh_parents_g <- c(
  paste( inputDir_g,   "12_Q18-W1579F-father-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_W1647F-father-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "14_W1719F-father-ICSI-fresh_Rep1.bismark.cov",   sep="/" ),
  paste( inputDir_g,   "12_Q19-W1579M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "13_Q23-W1647M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ),
  paste( inputDir_g,   "14_Q15-W1719M-mother-ICSI-fresh_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_fresh_parents_g <- rep( x="ICSI_fresh", times=length(Files_ICSI_fresh_parents_g) )
Tech_ICSI_fresh_parents_g <- mySampleID_ICSI_fresh_parents_g

for(i in c(1:length(mySampleID_ICSI_fresh_parents_g)) ) {
  mySampleID_ICSI_fresh_parents_g[i] = paste(mySampleID_ICSI_fresh_parents_g[i], i, sep="_")
}

myTreatment_ICSI_fresh_parents_g <- rep( x=3,  times=length(Files_ICSI_fresh_parents_g) )

Sex_ICSI_fresh_parents_g = rep( x="father",      times=length(Files_ICSI_fresh_parents_g) )  
for(i in c(1:length(Sex_ICSI_fresh_parents_g)) ) {
  Sex_ICSI_fresh_parents_g[i] = "father"
  if(i>=4) { Sex_ICSI_fresh_parents_g[i] = "mother" }
}



Files_ICSI_frozen_parents_g <- c(
  paste( inputDir_g,   "15_Q7-W871F-father-ICSI-frozen_Rep1.bismark.cov",     sep="/" ),
  paste( inputDir_g,   "15_Q10-W871M-mother-ICSI-frozen_Rep1.bismark.cov",    sep="/" ) 
)

mySampleID_ICSI_frozen_parents_g <- rep( x="ICSI_frozen", times=length(Files_ICSI_frozen_parents_g) )
Tech_ICSI_frozen_parents_g <- mySampleID_ICSI_frozen_parents_g

for(i in c(1:length(mySampleID_ICSI_frozen_parents_g)) ) {
  mySampleID_ICSI_frozen_parents_g[i] = paste(mySampleID_ICSI_frozen_parents_g[i], i, sep="_")
}

myTreatment_ICSI_frozen_parents_g <- rep( x=4,  times=length(Files_ICSI_frozen_parents_g) )

Sex_ICSI_frozen_parents_g = rep( x="father",      times=length(Files_ICSI_frozen_parents_g) )  
for(i in c(1:length(Sex_ICSI_frozen_parents_g)) ) {
  Sex_ICSI_frozen_parents_g[i] = "father"
  if(i>=2) { Sex_ICSI_frozen_parents_g[i] = "mother" }
}



Files_All_vector_parents_g <- c(
  Files_NC_parents_g,
  Files_IVF_fresh_parents_g,
  Files_IVF_frozen_parents_g,
  Files_ICSI_fresh_parents_g,
  Files_ICSI_frozen_parents_g  
)
Files_All_list_parents_g <- as.list( Files_All_vector_parents_g )


mySampleID_All_vector_parents_g <- c( 
  mySampleID_NC_parents_g,
  mySampleID_IVF_fresh_parents_g,
  mySampleID_IVF_frozen_parents_g,
  mySampleID_ICSI_fresh_parents_g,
  mySampleID_ICSI_frozen_parents_g  
)
mySampleID_All_list_parents_g <- as.list( mySampleID_All_vector_parents_g )


myTreatment_All_vector_parents_g <- c( 
  myTreatment_NC_parents_g,
  myTreatment_IVF_fresh_parents_g,
  myTreatment_IVF_frozen_parents_g,
  myTreatment_ICSI_fresh_parents_g,
  myTreatment_ICSI_frozen_parents_g 
)       
myTreatment_All_list_parents_g <- as.list( myTreatment_All_vector_parents_g )


mySex_All_vector_parents_g <- c( 
  Sex_NC_parents_g,
  Sex_IVF_fresh_parents_g,
  Sex_IVF_frozen_parents_g,
  Sex_ICSI_fresh_parents_g,
  Sex_ICSI_frozen_parents_g 
)       
mySex_All_list_parents_g <- as.list( mySex_All_vector_parents_g )


myTech_All_vector_parents_g <- c( 
  Tech_NC_parents_g,
  Tech_IVF_fresh_parents_g,
  Tech_IVF_frozen_parents_g,
  Tech_ICSI_fresh_parents_g,
  Tech_ICSI_frozen_parents_g 
)       
myTech_All_list_parents_g <- as.list( myTech_All_vector_parents_g )


## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple
mySex_All_shape_parents_g  = c( "boy"=16, "girl"=17, "father"=1, "mother"=2 ) 
myTech_All_color_parents_g = c( "NC"="black", "IVF_fresh"="blue", "IVF_frozen"="green",  "ICSI_fresh"="red", "ICSI_frozen"="purple" )
mySex_All_shape_parents_g
myTech_All_color_parents_g

MySex_Shape_parents_g <- function(  mySex_vector   ) {
  mySex_shape2  = mySex_vector
  for(i in c(1:length(mySex_shape2)) ) {
    if(mySex_shape2[i] == "boy")    { mySex_shape2[i] = c("boy"=16)   }
    if(mySex_shape2[i] == "girl")   { mySex_shape2[i] = c("girl"=17)  }
    if(mySex_shape2[i] == "father") { mySex_shape2[i] = c("father"=1) }
    if(mySex_shape2[i] == "mother") { mySex_shape2[i] = c("mother"=2) }
  }
  mySex_shape2 = as.numeric(mySex_shape2)
  names(mySex_shape2) = mySex_vector
  return(mySex_shape2)
}

MyTech_color_parents_g <- function(  myTech_vector  ) {
  myTech_color  = myTech_vector
  for(i in c(1:length(myTech_color)) ) {
    if(myTech_color[i] == "NC")          { myTech_color[i] = c("NC"="black") }
    if(myTech_color[i] == "IVF_fresh")   { myTech_color[i] = c("IVF_fresh"="blue") }
    if(myTech_color[i] == "IVF_frozen")  { myTech_color[i] = c("IVF_frozen"="green") }
    if(myTech_color[i] == "ICSI_fresh")  { myTech_color[i] = c("ICSI_fresh"="red") }
    if(myTech_color[i] == "ICSI_frozen") { myTech_color[i] = c("ICSI_frozen"="purple") }
  }
  names(myTech_color) = myTech_vector
  return(myTech_color)
}

MySex_Shape_parents_g( mySex_All_vector_parents_g )
MyTech_color_parents_g( myTech_All_vector_parents_g )
## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=blue,   IVF-frozen=green,  ICSI-fresh=red, ICSI-frozen=purple


length( Files_All_vector_parents_g )
length( Files_All_list_parents_g )
length( mySampleID_All_vector_parents_g )
length( mySampleID_All_list_parents_g )
length( myTreatment_All_vector_parents_g )
length( myTreatment_All_list_parents_g )
length( mySex_All_vector_parents_g )
length( mySex_All_list_parents_g )
length( myTech_All_vector_parents_g )
length( myTech_All_list_parents_g )
##################################################################################################################







##################################################################################################################
length( Files_All_vector_children_g )
length( Files_All_list_children_g )
length( mySampleID_All_vector_children_g )
length( mySampleID_All_list_children_g )
length( myTreatment_All_vector_children_g )
length( myTreatment_All_list_children_g )
length( mySex_All_vector_children_g )
length( mySex_All_list_children_g )
length( myTech_All_vector_children_g )
length( myTech_All_list_children_g )

length( Files_All_vector_parents_g )
length( Files_All_list_parents_g )
length( mySampleID_All_vector_parents_g )
length( mySampleID_All_list_parents_g )
length( myTreatment_All_vector_parents_g )
length( myTreatment_All_list_parents_g )
length( mySex_All_vector_parents_g )
length( mySex_All_list_parents_g )
length( myTech_All_vector_parents_g )
length( myTech_All_list_parents_g )



Files_All_vector_merge_g       = c(Files_All_vector_children_g,       Files_All_vector_parents_g)
mySampleID_All_vector_merge_g  = c(mySampleID_All_vector_children_g,  mySampleID_All_vector_parents_g)
myTreatment_All_vector_merge_g = c(myTreatment_All_vector_children_g, myTreatment_All_vector_parents_g)
mySex_All_vector_merge_g       = c(mySex_All_vector_children_g,       mySex_All_vector_parents_g)
myTech_All_vector_merge_g      = c(myTech_All_vector_children_g,      myTech_All_vector_parents_g)

Files_All_list_merge_g        =  as.list( Files_All_vector_merge_g ) 
mySampleID_All_list_merge_g   =  as.list( mySampleID_All_vector_merge_g ) 
myTreatment_All_list_merge_g  =  as.list( myTreatment_All_vector_merge_g ) 
mySex_All_list_merge_g        =  as.list( mySex_All_vector_merge_g )   
myTech_All_list_merge_g       =  as.list( myTech_All_vector_merge_g ) 



length( Files_All_vector_merge_g )
length( Files_All_list_merge_g )
length( mySampleID_All_vector_merge_g )
length( mySampleID_All_list_merge_g )
length( myTreatment_All_vector_merge_g )
length( myTreatment_All_list_merge_g )
length( mySex_All_vector_merge_g )
length( mySex_All_list_merge_g )
length( myTech_All_vector_merge_g )
length( myTech_All_list_merge_g )

##################################################################################################################







## Read the files about children
##################################################################################################################
myOutDir_sub1_children_g = paste(outDir_g, "/1A_ReadRawFiles_children",  sep="") 
if( ! file.exists(myOutDir_sub1_children_g) ) { dir.create(myOutDir_sub1_children_g, recursive = TRUE) }

sink( file=paste(myOutDir_sub1_children_g, "1_length-variables.txt", sep="/") )
length( Files_All_vector_children_g )
length( Files_All_list_children_g )
length( mySampleID_All_vector_children_g )
length( mySampleID_All_list_children_g )
length( myTreatment_All_vector_children_g )
length( myTreatment_All_list_children_g )
length( mySex_All_vector_children_g )
length( mySex_All_list_children_g )
length( myTech_All_vector_children_g )
length( myTech_All_list_children_g )
print( "#################### Files_All_vector_children_g ####################" )
print( Files_All_vector_children_g )
print( "#################### Files_All_list_children_g ####################" )
print( Files_All_list_children_g )
print( "#################### mySampleID_All_vector_children_g ####################" )
print( mySampleID_All_vector_children_g )
print( "#################### mySampleID_All_list_children_g ####################" )
print( mySampleID_All_list_children_g )
print( "#################### myTreatment_All_vector_children_g ####################" )
print( myTreatment_All_vector_children_g )
print( "#################### myTreatment_All_list_children_g ####################" )
print( myTreatment_All_list_children_g )
print( "#################### mySex_All_vector_children_g ####################" )
print( mySex_All_vector_children_g )
print( "#################### mySex_All_list_children_g ####################" )
print( mySex_All_list_children_g )
print( "#################### myTech_All_vector_children_g ####################" )
print( myTech_All_vector_children_g )
print( "#################### myTech_All_list_children_g ####################" )
print( myTech_All_list_children_g )
sink()


# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1_children_g, "2_theLog-of-read-inputFiles.txt", sep="/") )
myobj_children_g = methRead( Files_All_list_children_g,
                             sample.id = mySampleID_All_list_children_g,
                             assembly  = "hg38",
                             treatment = myTreatment_All_vector_children_g,
                             context   = "CpG",
                             pipeline  = "bismarkCoverage",
                             mincov    = 1,       ## >= n
                             header    = FALSE
)
sink()


sink( file=paste(myOutDir_sub1_children_g, "3_all-rawFiles.txt", sep="/") )
print(Files_All_vector_children_g)
print("#########################")
print(myobj_children_g)
sink()


sink( file=paste(myOutDir_sub1_children_g, "4_dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_children_g)) ) {
  print( "######################" )
  print(   Files_All_vector_children_g[i]  )
  print(   dim(myobj_children_g[[i]])  )
}
sink()


sink( file=paste(myOutDir_sub1_children_g, "5_dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_children_g)) ) {
  print(   dim(myobj_children_g[[i]])  )
}
sink()


continue_on_error_children_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_children_g) 


myobj_nor_children_g <- myobj_children_g
##myobj_nor_children_g <- normalizeCoverage(myobj_children_g)
##################################################################################################################





## Read the files about parents
##################################################################################################################
myOutDir_sub1_parents_g = paste(outDir_g, "/1B_ReadRawFiles_parents",  sep="") 
if( ! file.exists(myOutDir_sub1_parents_g) ) { dir.create(myOutDir_sub1_parents_g, recursive = TRUE) }

sink( file=paste(myOutDir_sub1_parents_g, "1_length-variables.txt", sep="/") )
length( Files_All_vector_parents_g )
length( Files_All_list_parents_g )
length( mySampleID_All_vector_parents_g )
length( mySampleID_All_list_parents_g )
length( myTreatment_All_vector_parents_g )
length( myTreatment_All_list_parents_g )
length( mySex_All_vector_parents_g )
length( mySex_All_list_parents_g )
length( myTech_All_vector_parents_g )
length( myTech_All_list_parents_g )
print( "#################### Files_All_vector_parents_g ####################" )
print( Files_All_vector_parents_g )
print( "#################### Files_All_list_parents_g ####################" )
print( Files_All_list_parents_g )
print( "#################### mySampleID_All_vector_parents_g ####################" )
print( mySampleID_All_vector_parents_g )
print( "#################### mySampleID_All_list_parents_g ####################" )
print( mySampleID_All_list_parents_g )
print( "#################### myTreatment_All_vector_parents_g ####################" )
print( myTreatment_All_vector_parents_g )
print( "#################### myTreatment_All_list_parents_g ####################" )
print( myTreatment_All_list_parents_g )
print( "#################### mySex_All_vector_parents_g ####################" )
print( mySex_All_vector_parents_g )
print( "#################### mySex_All_list_parents_g ####################" )
print( mySex_All_list_parents_g )
print( "#################### myTech_All_vector_parents_g ####################" )
print( myTech_All_vector_parents_g )
print( "#################### myTech_All_list_parents_g ####################" )
print( myTech_All_list_parents_g )
sink()


# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1_parents_g, "2_theLog-of-read-inputFiles.txt", sep="/") )
myobj_parents_g = methRead( Files_All_list_parents_g,
                            sample.id = mySampleID_All_list_parents_g,
                            assembly  = "hg38",
                            treatment = myTreatment_All_vector_parents_g,
                            context   = "CpG",
                            pipeline  = "bismarkCoverage",
                            mincov    = 1,       ## >= n
                            header    = FALSE
)
sink()


sink( file=paste(myOutDir_sub1_parents_g, "3_all-rawFiles.txt", sep="/") )
print(Files_All_vector_parents_g)
print("#########################")
print(myobj_parents_g)
sink()


sink( file=paste(myOutDir_sub1_parents_g, "4_dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_parents_g)) ) {
  print( "######################" )
  print(   Files_All_vector_parents_g[i]  )
  print(   dim(myobj_parents_g[[i]])  )
}
sink()


sink( file=paste(myOutDir_sub1_parents_g, "5_dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_parents_g)) ) {
  print(   dim(myobj_parents_g[[i]])  )
}
sink()


continue_on_error_parents_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_parents_g) 


myobj_nor_parents_g <- myobj_parents_g
##myobj_nor_parents_g <- normalizeCoverage(myobj_parents_g)
##################################################################################################################






## Read the files about parents
##################################################################################################################
myOutDir_sub1_merge_g = paste(outDir_g, "/1C_ReadRawFiles_allSamples",  sep="") 
if( ! file.exists(myOutDir_sub1_merge_g) ) { dir.create(myOutDir_sub1_merge_g, recursive = TRUE) }

sink( file=paste(myOutDir_sub1_merge_g, "1_length-variables.txt", sep="/") )
length( Files_All_vector_merge_g )
length( Files_All_list_merge_g )
length( mySampleID_All_vector_merge_g )
length( mySampleID_All_list_merge_g )
length( myTreatment_All_vector_merge_g )
length( myTreatment_All_list_merge_g )
length( mySex_All_vector_merge_g )
length( mySex_All_list_merge_g )
length( myTech_All_vector_merge_g )
length( myTech_All_list_merge_g )
print( "#################### Files_All_vector_merge_g ####################" )
print( Files_All_vector_merge_g )
print( "#################### Files_All_list_merge_g ####################" )
print( Files_All_list_merge_g )
print( "#################### mySampleID_All_vector_merge_g ####################" )
print( mySampleID_All_vector_merge_g )
print( "#################### mySampleID_All_list_merge_g ####################" )
print( mySampleID_All_list_merge_g )
print( "#################### myTreatment_All_vector_merge_g ####################" )
print( myTreatment_All_vector_merge_g )
print( "#################### myTreatment_All_list_merge_g ####################" )
print( myTreatment_All_list_merge_g )
print( "#################### mySex_All_vector_merge_g ####################" )
print( mySex_All_vector_merge_g )
print( "#################### mySex_All_list_merge_g ####################" )
print( mySex_All_list_merge_g )
print( "#################### myTech_All_vector_merge_g ####################" )
print( myTech_All_vector_merge_g )
print( "#################### myTech_All_list_merge_g ####################" )
print( myTech_All_list_merge_g )
sink()


# read the files to a methylRawList object: myobj
sink( file=paste(myOutDir_sub1_merge_g, "2_theLog-of-read-inputFiles.txt", sep="/") )
myobj_merge_g = methRead( Files_All_list_merge_g,
                          sample.id = mySampleID_All_list_merge_g,
                          assembly  = "hg38",
                          treatment = myTreatment_All_vector_merge_g,
                          context   = "CpG",
                          pipeline  = "bismarkCoverage",
                          mincov    = 1,       ## >= n
                          header    = FALSE
)
sink()


sink( file=paste(myOutDir_sub1_merge_g, "3_all-rawFiles.txt", sep="/") )
print(Files_All_vector_merge_g)
print("#########################")
print(myobj_merge_g)
sink()


sink( file=paste(myOutDir_sub1_merge_g, "4_dimensions-of-eachFile-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_merge_g)) ) {
  print( "######################" )
  print(   Files_All_vector_merge_g[i]  )
  print(   dim(myobj_merge_g[[i]])  )
}
sink()


sink( file=paste(myOutDir_sub1_merge_g, "5_dimensions-of-eachCov.txt", sep="/")  )
for( i in c(1:length(Files_All_vector_merge_g)) ) {
  print(   dim(myobj_merge_g[[i]])  )
}
sink()


continue_on_error_merge_g <- function()  {
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'")
}
# This is the key option
options(error=continue_on_error_merge_g) 


myobj_nor_merge_g <- myobj_merge_g
##myobj_nor_merge_g <- normalizeCoverage(myobj_merge_g)
##################################################################################################################







##################################################################################################################
myOutDir_sub2_g = paste(outDir_g, "/2_SD_CV_children_parents",  sep="") 
if( ! file.exists(myOutDir_sub2_g) ) { dir.create(myOutDir_sub2_g, recursive = TRUE) }

tiles_2two_children_g = tileMethylCounts( myobj_nor_children_g,   win.size=1000,   step.size=1000,   cov.bases = 3  )    
meth_2two_children_g  = unite( tiles_2two_children_g, destrand=FALSE, mc.cores=16   )   ## 100% overlap
mat_2two_children_g   = percMethylation( meth_2two_children_g )
head(mat_2two_children_g)
dim(mat_2two_children_g)

sd_2two_children_g = rowSds(mat_2two_children_g)
length(sd_2two_children_g)
head(sd_2two_children_g)

mean_2two_children_g = rowMeans2(mat_2two_children_g)
length(mean_2two_children_g)
head(mean_2two_children_g)

cv_2two_children_g = sd_2two_children_g/mean_2two_children_g
length(cv_2two_children_g)
head(cv_2two_children_g)

sd_mean_cv_2two_children_g  <- cbind(getData(meth_2two_children_g)[,1:4], sd_2two_children_g, mean_2two_children_g, cv_2two_children_g)               
head(sd_mean_cv_2two_children_g)


write.table(meth_2two_children_g , 
            file = paste(myOutDir_sub2_g,   "1A_meth-tiles_children.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
write.table(mat_2two_children_g , 
            file = paste(myOutDir_sub2_g,   "1B_mat-tiles_children.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
write.table(sd_mean_cv_2two_children_g , 
            file = paste(myOutDir_sub2_g,   "1C_SD-Mean-CV_children.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")





##############
tiles_2two_parents_g = tileMethylCounts( myobj_nor_parents_g,   win.size=1000,   step.size=1000,   cov.bases = 3  )    
meth_2two_parents_g  = unite( tiles_2two_parents_g, destrand=FALSE, mc.cores=16   )   ## 100% overlap
mat_2two_parents_g   = percMethylation( meth_2two_parents_g )
head(mat_2two_parents_g)
dim(mat_2two_parents_g)

sd_2two_parents_g = rowSds(mat_2two_parents_g)
length(sd_2two_parents_g)
head(sd_2two_parents_g)

mean_2two_parents_g = rowMeans2(mat_2two_parents_g)
length(mean_2two_parents_g)
head(mean_2two_parents_g)

cv_2two_parents_g = sd_2two_parents_g/mean_2two_parents_g
length(cv_2two_parents_g)
head(cv_2two_parents_g)


sd_mean_cv_2two_parents_g  <- cbind(getData(meth_2two_parents_g)[,1:4], sd_2two_parents_g, mean_2two_parents_g, cv_2two_parents_g)               
head(sd_mean_cv_2two_parents_g)


write.table(meth_2two_parents_g , 
            file = paste(myOutDir_sub2_g,   "2A_meth-tiles_parents.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
write.table(mat_2two_parents_g , 
            file = paste(myOutDir_sub2_g,   "2B_mat-tiles_parents.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
write.table(sd_mean_cv_2two_parents_g , 
            file = paste(myOutDir_sub2_g,   "2C_SD-Mean-CV_parents.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")



#################
dim(sd_mean_cv_2two_children_g)
dim(sd_mean_cv_2two_parents_g)

my_wdata_sd_g = data.frame( sd1 = c(sd_2two_children_g, sd_2two_parents_g),
                            type1 = factor( c(rep("children" , times=length(sd_2two_children_g)),
                                              rep("parents" ,  times=length(sd_2two_parents_g)) ) )
                          )

pdf( file=paste(myOutDir_sub2_g, "3A_StandardDeviation-tiles.pdf", sep="/") , width=3, height=3 )
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density")  
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 20)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 15)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 12.5)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 10)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 8)
dev.off()



pdf( file=paste(myOutDir_sub2_g, "3B_StandardDeviation-fill-tiles.pdf", sep="/") , width=3, height=3 )
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density")  
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 20)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 15)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 12.5)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 10)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 8)
dev.off()







pdf( file=paste(myOutDir_sub2_g, "3C_StandardDeviation-fill-alpha0.1-tiles.pdf", sep="/") , width=3, height=3 )
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density")  
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 20)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 15)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 12.5)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 10)
ggdensity(my_wdata_sd_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Standard Deviation", 
          xlab = "Standard Deviation (%)", ylab = "Density") + xlim(0, 8)
dev.off()



 

pdf( file=paste(myOutDir_sub2_g, "4A_StandardDeviation-CDF-tiles.pdf", sep="/"), width=3, height=3 )
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=1, 
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF")  
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 20)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 15)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 12)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 10)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 8)
dev.off()


 
pdf( file=paste(myOutDir_sub2_g, "4B_StandardDeviation-CDF-tiles.pdf", sep="/"), width=3, height=3 )
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2, 
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF")  
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 20)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 15)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 12)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 10)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 8)
dev.off()



pdf( file=paste(myOutDir_sub2_g, "4C_StandardDeviation-CDF-tiles.pdf", sep="/"), width=3, height=3 )
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",   
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF")  
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 20)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 15)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 12)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 10)
ggecdf(my_wdata_sd_g, x = "sd1", color = "type1", linetype = "type1",    
       title = "Distribution of Standard Deviation", 
       xlab = "Standard Deviation", ylab = "CDF") + xlim(0, 8)
dev.off()



png( file=paste(myOutDir_sub2_g, "5_StandardDeviation-QQplot-tiles.png", sep="/") )
ggqqplot(my_wdata_sd_g, x = "sd1", color = "type1")
dev.off()




##############
my_wdata_cv_g = data.frame( sd1 = c(cv_2two_children_g, cv_2two_parents_g),
                            type1 = factor( c(rep("children" , times=length(cv_2two_children_g)),
                                              rep("parents" ,  times=length(cv_2two_parents_g)) ) )
)


pdf( file=paste(myOutDir_sub2_g, "6A_CoefficientOfVariation-tiles.pdf", sep="/") , width=3, height=3 )
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density")  
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 1)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.5)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.3)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.2)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1", 
          alpha = 0.5, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.15)
dev.off()



pdf( file=paste(myOutDir_sub2_g, "6B_CoefficientOfVariation-fill-tiles.pdf", sep="/") , width=3, height=3 )
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density")  
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 1)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.5)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.3)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.2)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.2, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.15)
dev.off()







pdf( file=paste(myOutDir_sub2_g, "6C_CoefficientOfVariation-fill-alpha0.1-tiles.pdf", sep="/") , width=3, height=3 )
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density")  
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 1)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.5)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.3)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.2)
ggdensity(my_wdata_cv_g, x = "sd1",  add = "mean", color = "type1",  fill = "type1", 
          alpha = 0.1, title = "Distribution of Coefficient of Variation", 
          xlab = "Coefficient of Variation (%)", ylab = "Density") + xlim(0, 0.15)
dev.off()





pdf( file=paste(myOutDir_sub2_g, "7A_CoefficientOfVariation-CDF-tiles.pdf", sep="/"), width=3, height=3 )
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=1, 
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF")  
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 1)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.5)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.3)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.2)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=1,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.15)
dev.off()



pdf( file=paste(myOutDir_sub2_g, "7B_CoefficientOfVariation-CDF-tiles.pdf", sep="/"), width=3, height=3 )
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2, 
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF")  
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 1)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.5)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.3)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.2)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",  size=0.2,  
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.15)
dev.off()



pdf( file=paste(myOutDir_sub2_g, "7C_CoefficientOfVariation-CDF-tiles.pdf", sep="/"), width=3, height=3 )
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",   
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF")  
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 1)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.5)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.3)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",     
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.2)
ggecdf(my_wdata_cv_g, x = "sd1", color = "type1", linetype = "type1",    
       title = "Distribution of Coefficient of Variation", 
       xlab = "Coefficient of Variation", ylab = "CDF") + xlim(0, 0.15)
dev.off()



png( file=paste(myOutDir_sub2_g, "8_CoefficientOfVariation-QQplot-tiles.png", sep="/") )
ggqqplot(my_wdata_cv_g, x = "sd1", color = "type1")
dev.off()





##############
dim(sd_mean_cv_2two_children_g)
dim(sd_mean_cv_2two_parents_g)

myAnnotation_1( includeRegions2 = sd_mean_cv_2two_children_g, 
                path2 =  paste(myOutDir_sub2_g, "annotation_children", sep="/")
              ) 

myAnnotation_1( includeRegions2 = sd_mean_cv_2two_parents_g, 
                path2 =  paste(myOutDir_sub2_g, "annotation_parents", sep="/")
) 


##################################################################################################################
 








##################################################################################################################
myOutDir_sub3_g = paste(outDir_g, "/3_SD_CV_allSamples",  sep="") 
if( ! file.exists(myOutDir_sub3_g) ) { dir.create(myOutDir_sub3_g, recursive = TRUE) }

tiles_2two_merge_g = tileMethylCounts( myobj_nor_merge_g,   win.size=1000,   step.size=1000,   cov.bases = 3  )    
meth_2two_merge_g  = unite( tiles_2two_merge_g, destrand=FALSE, mc.cores=16   )   ## 100% overlap
mat_2two_merge_g   = percMethylation( meth_2two_merge_g )
head(mat_2two_merge_g)
dim(mat_2two_merge_g)

sd_2two_merge_g = rowSds(mat_2two_merge_g)
length(sd_2two_merge_g)
head(sd_2two_merge_g)

mean_2two_merge_g = rowMeans2(mat_2two_merge_g)
length(mean_2two_merge_g)
head(mean_2two_merge_g)

cv_2two_merge_g = sd_2two_merge_g/mean_2two_merge_g
length(cv_2two_merge_g)
head(cv_2two_merge_g)

sd_mean_cv_2two_merge_g  <- cbind(getData(meth_2two_merge_g)[,1:4], sd_2two_merge_g, mean_2two_merge_g, cv_2two_merge_g)               
head(sd_mean_cv_2two_merge_g)


write.table(meth_2two_merge_g , 
            file = paste(myOutDir_sub3_g,   "1A_meth-tiles_merge.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
write.table(mat_2two_merge_g , 
            file = paste(myOutDir_sub3_g,   "1B_mat-tiles_merge.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
write.table(sd_mean_cv_2two_merge_g , 
            file = paste(myOutDir_sub3_g,   "1C_SD-Mean-CV_merge.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")



##############
dim(sd_mean_cv_2two_merge_g)

myAnnotation_1( includeRegions2 = sd_mean_cv_2two_merge_g, 
                path2 =  paste(myOutDir_sub3_g, "annotation_merge", sep="/")
) 

 
###########
dim(mat_2two_merge_g)

sd_2two_mergeParents_g = rowSds(mat_2two_merge_g[,27:52])
length(sd_2two_mergeParents_g)
head(sd_2two_mergeParents_g)

mean_2two_mergeParents_g = rowMeans2(mat_2two_merge_g[,27:52])
length(mean_2two_mergeParents_g)
head(mean_2two_mergeParents_g)

cv_2two_mergeParents_g = sd_2two_merge_g/mean_2two_merge_g
length(cv_2two_mergeParents_g)
head(cv_2two_mergeParents_g)

sd_mean_cv_2two_mergeParents_g  <- cbind(getData(meth_2two_merge_g)[,1:4], sd_2two_mergeParents_g, 
                                  mean_2two_mergeParents_g, cv_2two_mergeParents_g)               
head(sd_mean_cv_2two_mergeParents_g)

write.table(sd_mean_cv_2two_mergeParents_g , 
            file = paste(myOutDir_sub3_g,   "2_SD-Mean-CV_mergeParents.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")


###########
dim(mat_2two_merge_g)

sd_2two_mergeChildren_g = rowSds(mat_2two_merge_g[,1:26])
length(sd_2two_mergeChildren_g)
head(sd_2two_mergeChildren_g)

mean_2two_mergeChildren_g = rowMeans2(mat_2two_merge_g[,1:26])
length(mean_2two_mergeChildren_g)
head(mean_2two_mergeChildren_g)

cv_2two_mergeChildren_g = sd_2two_merge_g/mean_2two_merge_g
length(cv_2two_mergeChildren_g)
head(cv_2two_mergeChildren_g)

sd_mean_cv_2two_mergeChildren_g  <- cbind(getData(meth_2two_merge_g)[,1:4], sd_2two_mergeChildren_g, 
                                          mean_2two_mergeChildren_g, cv_2two_mergeChildren_g)               
head(sd_mean_cv_2two_mergeChildren_g)

write.table(sd_mean_cv_2two_mergeChildren_g , 
            file = paste(myOutDir_sub3_g,   "3_SD-Mean-CV_mergeChildren.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")


############################



MyScatterDiagram_1 <- function(vector2X,  vector2Y,  path2,   fileName2,   xLab2,   yLab2,  title2,  height2=4,  width2=4,  yMin2=0, yMax2=2, 
                               xMin2=0, xMax2=2,  alpha2=0.5, diffThres2=1.5, colours2=c("red", "blue", "purple"))
  
  MyScatterDiagram_2 <- function(vector2X,  vector2Y,  path2,   fileName2,   xLab2,   yLab2,  title2,  
                                 height2=4,  width2=4,  yMin2=0, yMax2=2,  xMin2=0, xMax2=2 )
    
    
##################################################################################################################






myTopPercent_1 <- function(matrixTemp2, sdParentsTemp2, pathTemp2, percentTemp2, sdMeanCV2) {
  #matrixTemp2 =  mat_2two_merge_g
  #sdParentsTemp2 = sd_2two_mergeParents_g
  #pathTemp2 = paste(myOutDir_sub3_g,   "rmParetns",  sep="/")
  #percentTemp2 = 0.01
  #sdMeanCV2 = sd_mean_cv_2two_mergeParents_g
  #dim(matrixTemp2)
  #length(sdParentsTemp2)
  if( ! file.exists(pathTemp2) ) { dir.create(pathTemp2, recursive = TRUE) }
  
  sdParentsTemp2_sorted <- sort( sdParentsTemp2, decreasing = TRUE )
  threshMe <- sdParentsTemp2_sorted[length(sdParentsTemp2)*percentTemp2]
  myBoole1 <- (sdParentsTemp2 <= threshMe)  ## will be kept
  
  write.table(sdParentsTemp2[!myBoole1] , 
              file = paste(pathTemp2,   "1A_topSD_removed.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
  
  write.table(sdParentsTemp2[myBoole1] , 
              file = paste(pathTemp2,   "1B_SD_kept.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
  
  write.table(matrixTemp2[!myBoole1,] , 
              file = paste(pathTemp2,   "2A_topSD_removed_matrix.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
  
  write.table(matrixTemp2[myBoole1,] , 
              file = paste(pathTemp2,   "2B_SD_kept_matrix.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
  
  write.table(sdMeanCV2[!myBoole1,] , 
              file = paste(pathTemp2,   "3A_topSD_removed_sdMeanCV.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
  
  write.table(sdMeanCV2[myBoole1,] , 
              file = paste(pathTemp2,   "3B_SD_kept_sdMeanCV.txt",  sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = FALSE, col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")
  
  
  matrixTemp2_children <- matrixTemp2[,1:26]
  matrixTemp2_parents  <- matrixTemp2[,27:52]
  
  matrixTemp2_children_removed <- matrixTemp2_children[!myBoole1,]
  matrixTemp2_children_kept   <- matrixTemp2_children[myBoole1,]
  
  matrixTemp2_parents_removed <- matrixTemp2_parents[!myBoole1,]
  matrixTemp2_parents_kept   <- matrixTemp2_parents[myBoole1,]
  
  children_removed_dir = paste(pathTemp2, "HierarchicalClustering_children_removed", sep="/")
  if( ! file.exists(children_removed_dir) ) { dir.create(children_removed_dir, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=matrixTemp2_children_removed,   
                                 path_temp1=children_removed_dir,   dataFrame_temp1=NULL  ) 
  
  children_kept_dir = paste(pathTemp2, "HierarchicalClustering_children_kept", sep="/")
  if( ! file.exists(children_kept_dir) ) { dir.create(children_kept_dir, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=matrixTemp2_children_kept,   
                                 path_temp1=children_kept_dir,   dataFrame_temp1=NULL  ) 
  
  
  parents_removed_dir = paste(pathTemp2, "HierarchicalClustering_parents_removed", sep="/")
  if( ! file.exists(parents_removed_dir) ) { dir.create(parents_removed_dir, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=matrixTemp2_parents_removed,   
                                 path_temp1=parents_removed_dir,   dataFrame_temp1=NULL  ) 
  
  parents_kept_dir = paste(pathTemp2, "HierarchicalClustering_parents_kept", sep="/")
  if( ! file.exists(parents_kept_dir) ) { dir.create(parents_kept_dir, recursive = TRUE) }
  myHierarchicalClustering_1_g(  mat_3three=matrixTemp2_parents_kept,   
                                 path_temp1=parents_kept_dir,   dataFrame_temp1=NULL  ) 
  
}






##################################################################################################################
myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "4_rmParetns_SD_top0.001",  sep="/"), 
               percentTemp2 = 0.001, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)


myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "5_rmParetns_SD_top0.01",  sep="/"), 
               percentTemp2 = 0.01, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)


myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "6_rmParetns_SD_top0.02",  sep="/"), 
               percentTemp2 = 0.02, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)



myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "7_rmParetns_SD_top0.03",  sep="/"), 
               percentTemp2 = 0.03, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)


myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "8_rmParetns_SD_top0.04",  sep="/"), 
               percentTemp2 = 0.04, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)


myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "9_rmParetns_SD_top0.05",  sep="/"), 
               percentTemp2 = 0.05, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)



myTopPercent_1(matrixTemp2 = mat_2two_merge_g, sdParentsTemp2 = sd_2two_mergeParents_g, 
               pathTemp2 = paste(outDir_g,   "10_rmParetns_SD_top0.1",  sep="/"), 
               percentTemp2 = 0.1, sdMeanCV2 = sd_mean_cv_2two_mergeParents_g)



##################################################################################################################







