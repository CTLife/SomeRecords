## All reps must be overlapped. (100% overlap)
## example:  
## Rscript  RRBS2_singleAllPups_3.R     11B_splitXY/100A-rmXY/0_cov5reads     2011B_splitXY/100A-rmXY/0_cov5reads/RRBS2_singleAllPups_3/NC-vs-ICSIfresh          
         



args <- commandArgs(TRUE)
print("args: ")
print(args[1])   
print(args[2])     
print("#############")

inputDir = args[1];     ## the path of input file
outDir   = args[2];     ## the path of output file

# inputDir =  "11B_splitXY/100A-rmXY/5_cov50reads"
# outDir   =  "2011B_splitXY/100A-rmXY/5_cov50reads/RRBS2_singleAllPups_3/NC-vs-ICSIfresh"




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
library("reshape2") 
library("psych")
library("minerva")


myFileLists <- list(
paste(inputDir, "27_E16C-boy-NC-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "40_E12C-girl-NC_Rep1.bismark.cov", sep="/"),
paste(inputDir, "41_E13C-girl-NC_Rep1.bismark.cov", sep="/"),

paste(inputDir, "34_W200C-boy-ICSI-fresh_Rep1.bismark.cov", sep="/"),
paste(inputDir, "35_W766C-boy-ICSI-fresh-merge_Rep1.bismark.cov", sep="/"),
paste(inputDir, "36_W928C-boy-ICSI-fresh-merge_Rep1.bismark.cov",  sep="/") ,
paste(inputDir, "48_W808C-girl-ICSI-fresh_Rep1.bismark.cov", sep="/"),
paste(inputDir, "49_W924C-girl-ICSI-fresh_Rep1.bismark.cov", sep="/"),
paste(inputDir, "50_W934C-girl-ICSI-fresh_Rep1.bismark.cov",  sep="/")  

)


mySampleID <- list( 
"NC1", "NC2", "NC3",    
"ICSI-fresh1",  "ICSI-fresh2", "ICSI-fresh3",     
"ICSI-fresh4",  "ICSI-fresh5", "ICSI-fresh6"      
)


myTreatment <- c( 
0, 0, 0,     
1, 1,  1 , 1, 1,  1 
)       


## labels of all samples
myType1 <-  as.vector( unlist(mySampleID) )


## boy=16, girl=17,  father=1, mother=2
## NC=black,   IVF-fresh=red,   IVF-frozen=purple,  ICSI-fresh=blue, ICSI-frozen=green

myType2        = c( 
rep(x="boy",       times=1) ,
rep(x="girl",      times=2) ,
rep(x="boy",       times=3) ,
rep(x="girl",      times=3)  
)
myType2_shape  = c( 
rep( c("boy"=16),  times=1 ) , 
rep( c("girl"=17), times=2 ) , 
rep( c("boy"=16),  times=3 ) , 
rep( c("girl"=17), times=3 )  
) 
myType2_shape2 =  c( "boy"=16,  "girl"=17 )  


myType3        = c( rep(x="NC", times=3) , rep(x="ICSI-fresh", times=6) )
myType3_color  = c( rep( c("NC"="black"), times=3) , rep( c("ICSI-fresh"="blue"), times=6) )
myType3_color2 = c("NC"="black" ,  "ICSI-fresh"="blue")

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


## df contains two columns, the first column (cond_col=1) is sample type, the second column (val_col=2) is value. (must be).
whisk_1 <- function(df, cond_col=1, val_col=2) {  
  require(reshape2)
  condname <- names(df)[cond_col]  ## save the name of the first column.
  names(df)[cond_col] <- "cond" 
  names(df)[val_col]  <- "value"
  b   <- boxplot(value~cond, data=df, plot=FALSE)   
  df2 <- cbind(as.data.frame(b$stats), c("min","lq","m","uq","max"))
  names(df2) <- c(levels(df$cond), "pos")
  df2 <- melt(df2, id="pos", variable.name="cond")
  df2 <- dcast(df2, cond~pos)   
  names(df2)[1] <- condname 
  print(df2)
  df2
}


MyCluster_1 <- function( mymeth1 ,  path1, file1, width1, height1, sdThres1 ) {   
  pdf( file=paste(path1, file1, sep="/") , width=width1, height=height1  )
  
  clusterSamples(mymeth1, dist="correlation", method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="euclidean",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="maximum",     method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )
  clusterSamples(mymeth1, dist="manhattan",   method="ward",     sd.filter=TRUE,   sd.threshold=sdThres1,        filterByQuantile=TRUE,   plot=TRUE )

  dev.off()
  dev.off()
  dev.off()
  dev.off()
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


## T test and Wilcoxon test for large datasets.  
MyHypothesisTest_5 <- function(vector1, vector2, file1) {
  sink(file=file1)
  
  print("######################## Apply continuity correction in the normal approximation for the p-value. ###############################################")
  print("##################################################################################################################################")
  wilcoxTest_1 <- wilcox.test(x=vector1, y=vector2 )
  print( wilcoxTest_1  )
  cat( "\n\nExact p-value:", wilcoxTest_1$p.value, "\n\n\n\n\n" )
  
  print("######################## T-test, var.equal=FALSE. ###############################################")
  print("##################################################################################################################################")
  tTest_1 <- t.test(x=vector1, y=vector2 )
  print( tTest_1  )
  cat( "\n\nExact p-value:", tTest_1$p.value, "\n\n\n\n\n" )
  
  sink() 
  sink()
}



## compare  probability  density
MyHistogram_1 <- function(vector2, sampleType2,  colours2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4,  xMin2=0,  xMax2=1.5,   yMin2=0,  yMax2=10) {
  vector2[vector2>xMax2] <- xMax2
  vector2[vector2<xMin2] <- xMin2
  dataframeB  <- data.frame( xAxis = vector2,  sampleType=sampleType2 )
  
  
  FigureTemp1 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
      geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.3  ) +
      scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  ylim(yMin2, yMax2) + 
      scale_x_continuous(limits=c(xMin2, xMax2)  ) +  
      #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
      MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=0, shape=1)) ) }
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-density-limitY",      sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  
  FigureTemp2 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
      geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.3  ) +
      scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  
      scale_x_continuous(limits=c(xMin2, xMax2)  ) +  
      #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
      MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=0, shape=1)) ) }
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-density",             sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  
  FigureTemp4 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
      geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE, alpha=0.0  ) +
      scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  
      scale_x_continuous(limits=c(xMin2, xMax2)  ) +  
      #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
      MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "-density2",            sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  
  FigureTemp5 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
      geom_line(stat="density", alpha=1.0 ) +
      scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  ylim(yMin2, yMax2) + 
      scale_x_continuous(limits=c(xMin2, xMax2)  ) +  
      #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
      MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "-density3-limitY",     sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  
  FigureTemp6 <- { ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType, colour=sampleType) )   +  xlab(xLab2) + ylab("Probability density") +  ggtitle(title2)  +  
      geom_line(stat="density", alpha=1.0 ) +
      scale_colour_manual( values=colours2   ) + scale_fill_manual( values = colours2) +  
      scale_x_continuous(limits=c(xMin2, xMax2)  ) +  
      #geom_hline(yintercept=-0.01, lty=1, col="white", size=0.6) +
      MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   +  guides( colour = guide_legend(override.aes = list(size=1.5, shape=1)) ) }
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp6,  path1=path2, fileName1=paste(fileName2, "-density3",            sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
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



MyDistributionLine_1 <- function(x2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=4,  width2=7, yintercept2=c(-10, 0, 10)   ) {  
  myTempFrame_1 <- data.frame( xAxis = c(1:length(x2)),      yAxis = sort(x2) )
  
  FigureTemp1 <- ggplot(myTempFrame_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
    geom_line(size=0.5, colour="red") + geom_hline(yintercept =yintercept2[1], colour="blue") + 
    geom_hline(yintercept = yintercept2[2], colour="blue") + geom_hline(yintercept = yintercept2[3], colour="blue") + 
    MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-1",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp1 <- ggplot(myTempFrame_1,   aes(x=xAxis, y=yAxis )  )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2) + 
    geom_line(size=0.5, colour="red") + geom_hline(yintercept =yintercept2[1], colour="blue") + 
    geom_hline(yintercept = yintercept2[2], colour="blue") + geom_hline(yintercept = yintercept2[3], colour="blue") + 
    MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )    
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-2-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)        
}  



## x and y are discrete, x is factor, y is numeric.
MyHeatmap_1 <- function(matrix2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, midpoint2=0.5,  limits2=c(0, 1)  ) {  
  heatmap_1 <- melt(matrix2)
  
  FigureTemp1 <- ggplot(data=heatmap_1, aes(y=as.numeric(Var1), x=as.factor(Var2), fill=value) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
    geom_tile() + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = midpoint2,  limits=limits2 )  + 
    scale_x_discrete(expand = c(0, 0))  + scale_y_discrete(expand = c(0, 0)) +  #scale_y_reverse() +  coord_flip() + 
    MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "-1",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot(data=heatmap_1, aes(y=Var1, x=Var2, fill=value) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
    geom_tile() + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green", midpoint = midpoint2,  limits=limits2 )  + 
    #scale_x_discrete(expand = c(0, 0))  + scale_y_discrete(expand = c(0, 0)) +  #scale_y_reverse() +  coord_flip() + 
    MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-2",  sep="",  collapse=NULL),  height1=height2, width1=width2)
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
myOutDir_2two = paste(myOutDir, "/2-AllSites",  sep="");
if( ! file.exists(myOutDir_2two) ) { dir.create(myOutDir_2two, recursive = TRUE) }


tiles_2two = tileMethylCounts( myobj, win.size=100, step.size=100) ## 100 bin
meth_2two  = unite( tiles_2two, destrand=FALSE   )   ## 100% overlap
mat_2two   = percMethylation( meth_2two )




#####
myOutDir_2two_sub1 = paste(myOutDir_2two, "/1-stats-Files",  sep="");
if( ! file.exists(myOutDir_2two_sub1) ) { dir.create(myOutDir_2two_sub1, recursive = TRUE) }

sink( file=paste(myOutDir_2two_sub1 , "0-dimensions-2kbBin.txt", sep="/")  )
print( tiles_2two )
print("#########dimensions:")
print( dim(meth_2two)  )   
print( dim(mat_2two)   )
sink()


for( i in c(1:length(myFileLists)) ) {
  temp_file  = myFileLists[[i]]
  temp_file2 = strsplit(temp_file, split="/", fixed = TRUE, perl = FALSE) 
  temp_file3 = as.vector( unlist(temp_file2) )
  temp_file4 = temp_file3[ length(temp_file3) ]
  
  write.table(tiles_2two[[i]] , 
            file = paste(myOutDir_2two_sub1, "/1_",  temp_file4,  sep=""), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
}


write.table(meth_2two , 
            file = paste(myOutDir_2two_sub1,   "2A-meth-2kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")


write.table(mat_2two , 
            file = paste(myOutDir_2two_sub1,   "2B-mat-2kbBin.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = FALSE,  col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")







##################
myOutDir_2two_sub2 = paste(myOutDir_2two, "/2-oneClass",  sep="");
if( ! file.exists(myOutDir_2two_sub2) ) { dir.create(myOutDir_2two_sub2, recursive = TRUE) }

dim(mat_2two)
head(mat_2two)
vec_2two_values <-   array(mat_2two) 
vec_2two_sex    <-   rep(myType2, each = nrow(mat_2two), times = 1)
vec_2two_tech   <-   rep(myType3, each = nrow(mat_2two), times = 1)
vec_2two_name   <-   rep(colnames(mat_2two), each = nrow(mat_2two), times = 1)

length(vec_2two_values)
length(vec_2two_sex)
length(vec_2two_tech)
length(vec_2two_name)

vec_2two_values[1:100] 
vec_2two_sex[1:100] 
vec_2two_tech[1:100] 
vec_2two_name[1:100] 


sink(file = paste(myOutDir_2two_sub2, "1-violinPlot-byName.txt", sep="/"))
MyBoxViolinPlot_1(vector2=vec_2two_values,   sampleType2=vec_2two_name,  colours2="red",   
                  path2=myOutDir_2two_sub2,   fileName2="1-violinPlot-byName",    
                  title2="By name",  xLab2="All samples",  yLab2="Mehtylation level (%)",    
                  height2=3,   width2=1+length(colnames(mat_2two))/2,   Ymin2=0, Ymax2=100)
sink()  
  


sink(file = paste(myOutDir_2two_sub2, "2-violinPlot-bySex.txt", sep="/"))
MyBoxViolinPlot_1(vector2=vec_2two_values,   sampleType2=vec_2two_sex,  colours2="red",   
                  path2=myOutDir_2two_sub2,   fileName2="2-violinPlot-bySex",    
                  title2="By sex",  xLab2="All samples",  yLab2="Mehtylation level (%)",    
                  height2=3,   width2=2,   Ymin2=0, Ymax2=100)
sink()  
sex_calss1 <-  c(myType2=="boy")
sex_calss2 <-  c(myType2=="girl")
sex_calss1
sex_calss2
MyHypothesisTest_5( vector1=array(mat_2two[,sex_calss1]), vector2=array(mat_2two[,sex_calss2]), 
                    file1=paste(myOutDir_2two_sub2, "2-violinPlot-bySex.HypoTest.txt", sep="/") )

 

sink(file = paste(myOutDir_2two_sub2, "3-violinPlot-byTech.txt", sep="/"))
MyBoxViolinPlot_1(vector2=vec_2two_values,   sampleType2=vec_2two_tech,  colours2=myType3_color2,   
                  path2=myOutDir_2two_sub2,   fileName2="3-violinPlot-byTech",    
                  title2="By tech",  xLab2="All samples",  yLab2="Mehtylation level (%)",    
                  height2=3,   width2=2,   Ymin2=0, Ymax2=100)
sink()  
tech_calss1 <-  c(myTreatment==0)
tech_calss2 <-  c(myTreatment==1)
tech_calss1
tech_calss2
MyHypothesisTest_5( vector1=array(mat_2two[,tech_calss1]), vector2=array(mat_2two[,tech_calss2]), 
                    file1=paste(myOutDir_2two_sub2, "3-violinPlot-byTech.HypoTest.txt", sep="/") )


sink(file = paste(myOutDir_2two_sub2, "4-violinPlot-byTechSex.txt", sep="/"))
MyBoxViolinPlot_2(vector2=vec_2two_values,   sampleType2=vec_2two_sex,  
                  sampleType3=vec_2two_tech ,  colours2=myType3_color2,   
                  path2=myOutDir_2two_sub2,   fileName2="4-violinPlot-byTechSex",    
                  title2="By tech",  xLab2="All samples",  yLab2="Mehtylation level (%)",    
                  height2=3,   width2=4.5,   Ymin2=0, Ymax2=100)
sink()  
 


MyHistogram_1(vector2=vec_2two_values,   sampleType2=vec_2two_sex,  colours2=c("boy"="blue", "girl"="red"),   
                  path2=myOutDir_2two_sub2,   fileName2="5A-compareDensity-bySex",    
                  title2="By sex", xLab2="Mehtylation level (%)",    
                  height2=3,   width2=6,    xMin2=0,  xMax2=100,   yMin2=0,  yMax2=0.05)


MyHistogram_1(vector2=vec_2two_values,   sampleType2=vec_2two_tech,  colours2=myType3_color2,   
              path2=myOutDir_2two_sub2,   fileName2="5B-compareDensity-byTech",    
              title2="By tech", xLab2="Mehtylation level (%)",    
              height2=3,   width2=6,    xMin2=0,  xMax2=100,   yMin2=0,  yMax2=0.05)



MyScatterDiagram_1(vector2X = rowMeans(mat_2two[,sex_calss1]) ,   
              vector2Y = rowMeans(mat_2two[,sex_calss2]) ,
              path2=myOutDir_2two_sub2,   fileName2="6A-scaterPlot-bySex",    
              title2="By sex", xLab2="Mehtylation level (%, boy)",   yLab2="Mehtylation level (%, girl)",  
              height2=3,   width2=4,    xMin2=0,  xMax2=100,   yMin2=0,  yMax2=100,  alpha2=0.5, 
              diffThres2=10, colours2=c("red", "blue", "purple"))

MyScatterDiagram_1(vector2X = rowMeans(mat_2two[,tech_calss1]) ,   
                   vector2Y = rowMeans(mat_2two[,tech_calss2]) ,
                   path2=myOutDir_2two_sub2,   fileName2="6B-scaterPlot-byTech",    
                   title2="By tech", xLab2="Mehtylation level (%, NC)",   yLab2="Mehtylation level (%)",  
                   height2=3,   width2=4,    xMin2=0,  xMax2=100,   yMin2=0,  yMax2=100,  alpha2=0.5, 
                   diffThres2=10, colours2=c("red", "blue", "purple"))





MyScatterDiagram_2(vector2X = rowMeans(mat_2two[,sex_calss1]) ,   
                   vector2Y = rowMeans(mat_2two[,sex_calss2]) ,
                   path2=myOutDir_2two_sub2,   fileName2="7A-densityPlot-bySex",    
                   title2="By sex", xLab2="Mehtylation level (%, boy)",   yLab2="Mehtylation level (%, girl)",  
                   height2=3,   width2=4,    xMin2=0,  xMax2=100,   yMin2=0,  yMax2=100 
                   )

MyScatterDiagram_2(vector2X = rowMeans(mat_2two[,tech_calss1]) ,   
                   vector2Y = rowMeans(mat_2two[,tech_calss2]) ,
                   path2=myOutDir_2two_sub2,   fileName2="7B-densityPlot-byTech",    
                   title2="By tech", xLab2="Mehtylation level (%, NC)",   yLab2="Mehtylation level (%)",  
                   height2=3,   width2=4,    xMin2=0,  xMax2=100,   yMin2=0,  yMax2=100 
                   )







MyDistributionLine_1(x2 = rowMeans(mat_2two[,sex_calss1]) - rowMeans(mat_2two[,sex_calss2]) ,   
                   path2=myOutDir_2two_sub2,   fileName2="8A-diffLine-bySex",    
                   title2="By sex", xLab2="sites",   yLab2="difference (%)",   
                   height2=3,   width2=5,    Ymin2= -50,   Ymax2=50,  yintercept2=c(-10, 0, 10)
)

MyDistributionLine_1(x2 = rowMeans(mat_2two[,tech_calss1]) - rowMeans(mat_2two[,tech_calss2]),   
                   path2=myOutDir_2two_sub2,   fileName2="8B-diffLine-byTech",    
                   title2="By tech", xLab2="sites",   yLab2="difference (%)",  
                   height2=3,   width2=5,    Ymin2= -50,   Ymax2=50,  yintercept2=c(-10, 0, 10) 
)




index_2two_1 = order( rowMeans(mat_2two[,tech_calss1]) )
mat_2two_sorted = mat_2two[index_2two_1, ]
rownames(mat_2two_sorted) = c(1:length(index_2two_1))

MyHeatmap_1(matrix2= mat_2two_sorted,  path2=myOutDir_2two_sub2,  fileName2="9A-heatmap-1", 
            title2="heatmap",  xLab2="samples",  yLab2="sites",   Ymin2=0,   Ymax2=100,   height2=3,  width2=4, midpoint2=50,  limits2=c(0, 100)  )


 


#####
##################
myOutDir_2two_sub3 = paste(myOutDir_2two, "/3-threeClass",  sep="");
if( ! file.exists(myOutDir_2two_sub3) ) { dir.create(myOutDir_2two_sub3, recursive = TRUE) }


mat_2two_sub3       <-  mat_2two 
myLow_2two_sub3     <-  as.matrix( mat_2two_sub3<20 )  
myMedium_2two_sub3  <-  as.matrix( (mat_2two_sub3>=20) & (mat_2two_sub3<=80) )
myHigh_2two_sub3    <-  as.matrix( mat_2two_sub3>80 )
 
mat_2two_sub3[myLow_2two_sub3]    <- "1Low"
mat_2two_sub3[myMedium_2two_sub3] <- "2Medium"
mat_2two_sub3[myHigh_2two_sub3]   <- "3High"

 
dim(mat_2two)
head(mat_2two)
dim(mat_2two_sub3)
head(mat_2two_sub3)
vec_2two_sub3 <-   array(mat_2two_sub3) 

length(vec_2two_sub3) 
length(vec_2two_values)
length(vec_2two_sex)
length(vec_2two_tech)
length(vec_2two_name)

vec_2two_sub3[1:100]
vec_2two_values[1:100] 
vec_2two_sex[1:100] 
vec_2two_tech[1:100] 
vec_2two_name[1:100] 






sink(file = paste(myOutDir_2two_sub2, "5-violinPlot-byStrength.txt", sep="/"))
MyBoxViolinPlot_2(vector2=vec_2two_values,   sampleType2=vec_2two_strength,  
                  sampleType3=vec_2two_tech ,  colours2=myType3_color2,   
                  path2=myOutDir_2two_sub2,   fileName2="5-violinPlot-byStrength",    
                  title2="By tech",  xLab2="All samples",  yLab2="Mehtylation level (%)",   
                  height2=3,   width2=5.5,   Ymin2=0, Ymax2=100)
sink()  



sink(file = paste(myOutDir_2two_sub2, "5A-violinPlot-byStrength-Low.txt", sep="/"))
MyBoxViolinPlot_1(vector2=vec_2two_values[myLow_2two_a],   sampleType2=vec_2two_tech[myLow_2two_a] ,  colours2=myType3_color2,   
                  path2=myOutDir_2two_sub2,   fileName2="5A-violinPlot-byStrength-Low",    
                  title2="By tech",  xLab2="All samples",  yLab2="Mehtylation level (%)",   
                  height2=3,   width2=2,   Ymin2=0, Ymax2=20)
sink()  
 




DataFrame_2two_sub2 <- data.frame( yAxis=vec_2two_values, mySex=vec_2two_sex, myTech=vec_2two_tech, myName=vec_2two_name ) 




######################################################################################################################################################
######################################################################################################################################################


















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


































