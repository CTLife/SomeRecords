## example:  Rscript   Cluster-PCA.R     1-ReadRawFiles  10B_mathylationLevel_5reads_mergeOverlap.txt    10-5reads



args <- commandArgs(TRUE)
print("args: ")
print(args[1])   
print(args[2])   
print(args[3])   
print("#############")

inputDir   = args[1];
inputFile  = args[2];     ## input file
outDir     = args[2];     ## output path
#  inputDir   = "1-ReadRawFiles";
#  inputFile  = "10B_mathylationLevel_5reads_mergeOverlap.txt";     ## input file
#  outDir = "10-5reads"

if( ! file.exists(outDir) ) { dir.create(outDir, recursive = TRUE)  }



#####################################################################################################################
## reduce columns of a matrix by average the nearest columns.
reduceMatrixCol <- function(matrix_1, colNum_1=10 ) {     ## colNum_1: number of columns after reduced.
  allCol  <- ncol(matrix_1)
  allRow  <- nrow(matrix_1)
  matrix2 <- matrix(nrow = allRow, ncol = colNum_1)
  aveCol  <- allCol/colNum_1
  for(i  in  c(0:(colNum_1-1)) ) {
    start <- floor(aveCol*i+1)
    end   <- floor(aveCol*(i+1))
    index2 <- c( start:end )
    cat(index2, "\n")
    matrix2[,i+1] <- rowMeans(matrix_1[, index2])
  }
  return(matrix2)
}

## reduce rows of a matrix by average the nearest rows.   
reduceMatrixRow <- function(matrix_1, rowNum_1=10 ) {     ## rowNum_1: number of rows after reduced.
  allCol  <- ncol(matrix_1)
  allRow  <- nrow(matrix_1)
  matrix2 <- matrix(nrow = rowNum_1, ncol = allCol)
  aveRow  <- allRow/rowNum_1
  for(i  in  c(0:(rowNum_1-1)) )  {
    start <- floor(aveRow*i+1)
    end   <- floor(aveRow*(i+1))
    index2 <- c( start:end )
    cat(index2, "\n")
    matrix2[i+1,] <- colMeans(matrix_1[index2, ])
  }
  return(matrix2)
}



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



#####################################################################################################################






###########################################################################################
matrix_1 <- read.table(paste(inputDir, inputFile, sep="/") , header=TRUE, sep="\t",  quote = "",  comment.char = "") 
dim(matrix_1)
head(matrix_1)
colnames(matrix_1)

vector_1 <- colMeans(matrix_1)
length(vector_1)


myType2 <- c(
  "girl",
  "girl",
  
  "girl",
  "girl",
  
  "girl",
  "girl",
  
  "girl",
  "girl",
  
  "girl",
  "girl",
  
  "girl",
  "girl",
  
  "boy",
  "boy",
  
  "boy",
  "boy",
  
  "boy",
  "boy",
  
  "boy",
  "boy",
  
  "boy",
  "boy",
  
  "boy",
  "boy",
  
  "boy",
  "girl",
  
  "girl",
  "boy",
  
  "girl",
  "boy",
  
  "girl",
  "boy",
  
  "boy",
  "girl",
  
  "girl",
  "boy",
  
  "girl",
  "boy"
  
)

myType3 <- c(
  "NC",
  "NC",
  
  "NC",
  "NC",
  
  "NC",
  "NC",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "NC",
  "NC",
  
  "NC",
  "NC",
  
  "NC",
  "NC",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "NC",
  "NC",
  
  "NC",
  "NC",
  
  "NC",
  "NC",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh",
  
  "IVF-fresh",
  "IVF-fresh"
  
)

## boy=16, girl=17,  father=1, mother=2
## NC=cyan,   IVF-fresh=red, IVF-frozen=purple,  ICSI-fresh=blue, ICSI-frozen=green

myType2_shape <- c(
  "girl"=17,
  "girl"=17,
  
  "girl"=17,
  "girl"=17,
  
  "girl"=17,
  "girl"=17,
  
  "girl"=17,
  "girl"=17,
  
  "girl"=17,
  "girl"=17,
  
  "girl"=17,
  "girl"=17,
  
  "boy"=16,
  "boy"=16,
  
  "boy"=16,
  "boy"=16,
  
  "boy"=16,
  "boy"=16,
  
  "boy"=16,
  "boy"=16,
  
  "boy"=16,
  "boy"=16,
  
  "boy"=16,
  "boy"=16,
  
  "boy"=16,
  "girl"=17,
  
  "girl"=17,
  "boy"=16,
  
  "girl"=17,
  "boy"=16,
  
  "girl"=17,
  "boy"=16,
  
  "boy"=16,
  "girl"=17,
  
  "girl"=17,
  "boy"=16,
  
  "girl"=17,
  "boy"=16
  
)


myType3_color <- c(
  "NC"="cyan",
  "NC"="cyan",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "NC"="cyan",
  "NC"="cyan",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red",
  
  "IVF-fresh"="red",
  "IVF-fresh"="red"
  
)


print( length(myType2) )
print( length(myType3) )
print( length(myType2_shape) )
print( length(myType3_color) )


myDataFrame <- data.frame( yAxis=vector_1,  mySex=myType2, myTech=myType3 ) 

library(ggplot2)

FigureTemp2B <- ggplot(myDataFrame, aes(x=mySex , y=yAxis,  fill=myTech, colour=myTech ) ) +  
  geom_boxplot( position=position_dodge(0.9), outlier.size=0 , alpha=0  ) +  
  geom_jitter( aes(colour=myTech ), size=2,  alpha=0.6, position = position_jitterdodge(jitter.width=0.7, dodge.width=0.9 )  ) +
  stat_summary( position=position_dodge(0.9), alpha=0.6, 
                fun.y=mean,  color="black",  geom="point",  shape=19, size=2, show.legend = FALSE  ) + 
  xlab( "sex and age" ) + ylab( "CpG sites mehtylation level (%)" ) +  
  scale_fill_manual( values = c("NC"="cyan",   "IVF-fresh"="red", "IVF-frozen"="purple",  "ICSI-fresh"="blue", "ICSI-frozen"="green") ) +
  ggtitle( ">= 1 reads" )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 )  
ggsave( filename = paste(outDir,  "/", "1-boxPlot-sex-tech",  ".svg",  sep="",  collapse=NULL),   height=4, width=6,  dpi = 1200 )




########################################
dim(matrix_1)
matrix_bin_100bp   <- reduceMatrixRow(matrix_1, rowNum_1=30297 )
matrix_bin_500bp   <- reduceMatrixRow(matrix_bin_100bp, rowNum_1=6059 )
matrix_bin_1000bp  <- reduceMatrixRow(matrix_bin_500bp, rowNum_1=3029 )
matrix_bin_5000bp  <- reduceMatrixRow(matrix_bin_1000bp, rowNum_1=605 )
matrix_bin_10000bp <- reduceMatrixRow(matrix_bin_5000bp, rowNum_1=302 )
matrix_bin_50000bp <- reduceMatrixRow(matrix_bin_10000bp, rowNum_1=60 )

dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)




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
library(ggrepel)
library(scatterplot3d)




######################################### 
myOutDir_sub2 = paste(outDir, "/2-bin-100bp",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

myPCA2 <- prcomp( t(matrix_bin_100bp) )
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




myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = colnames(matrix_1)
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


#####################



dim(matrix_1)
dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)



library("cluster")
library("factoextra")
library("magrittr")


colnames(matrix_1)
colnames(matrix_bin_100bp) <- colnames(matrix_1)
colnames(matrix_bin_500bp) <- colnames(matrix_1)
colnames(matrix_bin_1000bp) <- colnames(matrix_1)
colnames(matrix_bin_5000bp) <- colnames(matrix_1)
colnames(matrix_bin_10000bp) <- colnames(matrix_1)
colnames(matrix_bin_50000bp) <- colnames(matrix_1)


res.dist1 <- get_dist( t(matrix_bin_100bp) , stand = TRUE, method = "euclidean")
res.dist2 <- get_dist( t(matrix_bin_100bp) , stand = TRUE, method = "manhattan")
res.dist3 <- get_dist( t(matrix_bin_100bp) , stand = TRUE, method = "spearman")
res.dist4 <- get_dist( t(matrix_bin_100bp) , stand = TRUE, method = "pearson")

pdf( file = paste(myOutDir_sub2, "10_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


pdf( file = paste(myOutDir_sub2, "11_Determining-optimal-number-of-clusters.pdf",  sep="/") )
fviz_nbclust(t(matrix_bin_100bp), hcut,   method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
fviz_nbclust(t(matrix_bin_100bp), kmeans, method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
dev.off() 


km.res1 <- kmeans( t(matrix_bin_100bp), 6, nstart = 25)
km.res2 <- kmeans( t(matrix_bin_100bp), 7, nstart = 25)

pdf( file = paste(myOutDir_sub2, "12_kmeans.pdf",  sep="/") )
fviz_cluster(km.res1, data =  t(matrix_bin_100bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
fviz_cluster(km.res2, data =  t(matrix_bin_100bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
dev.off() 


## the k-medoids/pam clustering 
pam.res1 <- pam(t(matrix_bin_100bp), 6)
pam.res2 <- pam(t(matrix_bin_100bp), 7)
pdf( file = paste(myOutDir_sub2, "13_k-medoids.pdf",  sep="/") )
fviz_cluster(pam.res1)
fviz_cluster(pam.res2)
dev.off() 





# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "complete" )   
res.hc2 <- hclust(res.dist2, method = "complete" )   
res.hc3 <- hclust(res.dist3, method = "complete" )   
res.hc4 <- hclust(res.dist4, method = "complete" )   
pdf( file = paste(myOutDir_sub2, "14A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 



# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D" )   
res.hc2 <- hclust(res.dist2, method = "ward.D" )   
res.hc3 <- hclust(res.dist3, method = "ward.D" )   
res.hc4 <- hclust(res.dist4, method = "ward.D" )   
pdf( file = paste(myOutDir_sub2, "14B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D2" )   
res.hc2 <- hclust(res.dist2, method = "ward.D2" )   
res.hc3 <- hclust(res.dist3, method = "ward.D2" )   
res.hc4 <- hclust(res.dist4, method = "ward.D2" )   
pdf( file = paste(myOutDir_sub2, "14C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




## Fuzzy clustering
res.fanny1 <- fanny( t(matrix_bin_100bp), 4)
res.fanny2 <- fanny( t(matrix_bin_100bp), 5)
res.fanny3 <- fanny( t(matrix_bin_100bp), 6)
res.fanny4 <- fanny( t(matrix_bin_100bp), 7)

pdf( file = paste(myOutDir_sub2, "15_Fuzzy clustering.pdf",  sep="/") )
fviz_cluster(res.fanny1, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny2, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny3, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny4, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
dev.off() 





















######################################### 
myOutDir_sub2 = paste(outDir, "/3-bin-500bp",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

myPCA2 <- prcomp( t(matrix_bin_500bp) )
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




myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = colnames(matrix_1)
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


#####################



dim(matrix_1)
dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)



library("cluster")
library("factoextra")
library("magrittr")


colnames(matrix_1)
colnames(matrix_bin_100bp) <- colnames(matrix_1)
colnames(matrix_bin_500bp) <- colnames(matrix_1)
colnames(matrix_bin_1000bp) <- colnames(matrix_1)
colnames(matrix_bin_5000bp) <- colnames(matrix_1)
colnames(matrix_bin_10000bp) <- colnames(matrix_1)
colnames(matrix_bin_50000bp) <- colnames(matrix_1)


res.dist1 <- get_dist( t(matrix_bin_500bp) , stand = TRUE, method = "euclidean")
res.dist2 <- get_dist( t(matrix_bin_500bp) , stand = TRUE, method = "manhattan")
res.dist3 <- get_dist( t(matrix_bin_500bp) , stand = TRUE, method = "spearman")
res.dist4 <- get_dist( t(matrix_bin_500bp) , stand = TRUE, method = "pearson")

pdf( file = paste(myOutDir_sub2, "10_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


pdf( file = paste(myOutDir_sub2, "11_Determining-optimal-number-of-clusters.pdf",  sep="/") )
fviz_nbclust(t(matrix_bin_500bp), hcut,   method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
fviz_nbclust(t(matrix_bin_500bp), kmeans, method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
dev.off() 


km.res1 <- kmeans( t(matrix_bin_500bp), 6, nstart = 25)
km.res2 <- kmeans( t(matrix_bin_500bp), 7, nstart = 25)

pdf( file = paste(myOutDir_sub2, "12_kmeans.pdf",  sep="/") )
fviz_cluster(km.res1, data =  t(matrix_bin_500bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
fviz_cluster(km.res2, data =  t(matrix_bin_500bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
dev.off() 


## the k-medoids/pam clustering 
pam.res1 <- pam(t(matrix_bin_500bp), 6)
pam.res2 <- pam(t(matrix_bin_500bp), 7)
pdf( file = paste(myOutDir_sub2, "13_k-medoids.pdf",  sep="/") )
fviz_cluster(pam.res1)
fviz_cluster(pam.res2)
dev.off() 





# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "complete" )   
res.hc2 <- hclust(res.dist2, method = "complete" )   
res.hc3 <- hclust(res.dist3, method = "complete" )   
res.hc4 <- hclust(res.dist4, method = "complete" )   
pdf( file = paste(myOutDir_sub2, "14A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 



# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D" )   
res.hc2 <- hclust(res.dist2, method = "ward.D" )   
res.hc3 <- hclust(res.dist3, method = "ward.D" )   
res.hc4 <- hclust(res.dist4, method = "ward.D" )   
pdf( file = paste(myOutDir_sub2, "14B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D2" )   
res.hc2 <- hclust(res.dist2, method = "ward.D2" )   
res.hc3 <- hclust(res.dist3, method = "ward.D2" )   
res.hc4 <- hclust(res.dist4, method = "ward.D2" )   
pdf( file = paste(myOutDir_sub2, "14C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




## Fuzzy clustering
res.fanny1 <- fanny( t(matrix_bin_500bp), 4)
res.fanny2 <- fanny( t(matrix_bin_500bp), 5)
res.fanny3 <- fanny( t(matrix_bin_500bp), 6)
res.fanny4 <- fanny( t(matrix_bin_500bp), 7)

pdf( file = paste(myOutDir_sub2, "15_Fuzzy clustering.pdf",  sep="/") )
fviz_cluster(res.fanny1, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny2, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny3, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny4, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
dev.off() 
























######################################### 
myOutDir_sub2 = paste(outDir, "/4-bin-1000bp",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

myPCA2 <- prcomp( t(matrix_bin_1000bp) )
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




myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = colnames(matrix_1)
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




#####################



dim(matrix_1)
dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)



library("cluster")
library("factoextra")
library("magrittr")


colnames(matrix_1)
colnames(matrix_bin_100bp) <- colnames(matrix_1)
colnames(matrix_bin_500bp) <- colnames(matrix_1)
colnames(matrix_bin_1000bp) <- colnames(matrix_1)
colnames(matrix_bin_5000bp) <- colnames(matrix_1)
colnames(matrix_bin_10000bp) <- colnames(matrix_1)
colnames(matrix_bin_50000bp) <- colnames(matrix_1)


res.dist1 <- get_dist( t(matrix_bin_1000bp) , stand = TRUE, method = "euclidean")
res.dist2 <- get_dist( t(matrix_bin_1000bp) , stand = TRUE, method = "manhattan")
res.dist3 <- get_dist( t(matrix_bin_1000bp) , stand = TRUE, method = "spearman")
res.dist4 <- get_dist( t(matrix_bin_1000bp) , stand = TRUE, method = "pearson")

pdf( file = paste(myOutDir_sub2, "10_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


pdf( file = paste(myOutDir_sub2, "11_Determining-optimal-number-of-clusters.pdf",  sep="/") )
fviz_nbclust(t(matrix_bin_1000bp), hcut,   method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
fviz_nbclust(t(matrix_bin_1000bp), kmeans, method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
dev.off() 


km.res1 <- kmeans( t(matrix_bin_1000bp), 6, nstart = 25)
km.res2 <- kmeans( t(matrix_bin_1000bp), 7, nstart = 25)

pdf( file = paste(myOutDir_sub2, "12_kmeans.pdf",  sep="/") )
fviz_cluster(km.res1, data =  t(matrix_bin_1000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
fviz_cluster(km.res2, data =  t(matrix_bin_1000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
dev.off() 


## the k-medoids/pam clustering 
pam.res1 <- pam(t(matrix_bin_1000bp), 6)
pam.res2 <- pam(t(matrix_bin_1000bp), 7)
pdf( file = paste(myOutDir_sub2, "13_k-medoids.pdf",  sep="/") )
fviz_cluster(pam.res1)
fviz_cluster(pam.res2)
dev.off() 





# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "complete" )   
res.hc2 <- hclust(res.dist2, method = "complete" )   
res.hc3 <- hclust(res.dist3, method = "complete" )   
res.hc4 <- hclust(res.dist4, method = "complete" )   
pdf( file = paste(myOutDir_sub2, "14A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 



# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D" )   
res.hc2 <- hclust(res.dist2, method = "ward.D" )   
res.hc3 <- hclust(res.dist3, method = "ward.D" )   
res.hc4 <- hclust(res.dist4, method = "ward.D" )   
pdf( file = paste(myOutDir_sub2, "14B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D2" )   
res.hc2 <- hclust(res.dist2, method = "ward.D2" )   
res.hc3 <- hclust(res.dist3, method = "ward.D2" )   
res.hc4 <- hclust(res.dist4, method = "ward.D2" )   
pdf( file = paste(myOutDir_sub2, "14C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




## Fuzzy clustering
res.fanny1 <- fanny( t(matrix_bin_1000bp), 4)
res.fanny2 <- fanny( t(matrix_bin_1000bp), 5)
res.fanny3 <- fanny( t(matrix_bin_1000bp), 6)
res.fanny4 <- fanny( t(matrix_bin_1000bp), 7)

pdf( file = paste(myOutDir_sub2, "15_Fuzzy clustering.pdf",  sep="/") )
fviz_cluster(res.fanny1, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny2, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny3, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny4, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
dev.off() 


































######################################### 
myOutDir_sub2 = paste(outDir, "/5-bin-5000bp",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

myPCA2 <- prcomp( t(matrix_bin_5000bp) )
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




myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = colnames(matrix_1)
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




#####################



dim(matrix_1)
dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)



library("cluster")
library("factoextra")
library("magrittr")


colnames(matrix_1)
colnames(matrix_bin_100bp) <- colnames(matrix_1)
colnames(matrix_bin_500bp) <- colnames(matrix_1)
colnames(matrix_bin_1000bp) <- colnames(matrix_1)
colnames(matrix_bin_5000bp) <- colnames(matrix_1)
colnames(matrix_bin_10000bp) <- colnames(matrix_1)
colnames(matrix_bin_50000bp) <- colnames(matrix_1)


res.dist1 <- get_dist( t(matrix_bin_5000bp) , stand = TRUE, method = "euclidean")
res.dist2 <- get_dist( t(matrix_bin_5000bp) , stand = TRUE, method = "manhattan")
res.dist3 <- get_dist( t(matrix_bin_5000bp) , stand = TRUE, method = "spearman")
res.dist4 <- get_dist( t(matrix_bin_5000bp) , stand = TRUE, method = "pearson")

pdf( file = paste(myOutDir_sub2, "10_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


pdf( file = paste(myOutDir_sub2, "11_Determining-optimal-number-of-clusters.pdf",  sep="/") )
fviz_nbclust(t(matrix_bin_5000bp), hcut,   method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
fviz_nbclust(t(matrix_bin_5000bp), kmeans, method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
dev.off() 


km.res1 <- kmeans( t(matrix_bin_5000bp), 6, nstart = 25)
km.res2 <- kmeans( t(matrix_bin_5000bp), 7, nstart = 25)

pdf( file = paste(myOutDir_sub2, "12_kmeans.pdf",  sep="/") )
fviz_cluster(km.res1, data =  t(matrix_bin_5000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
fviz_cluster(km.res2, data =  t(matrix_bin_5000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
dev.off() 


## the k-medoids/pam clustering 
pam.res1 <- pam(t(matrix_bin_5000bp), 6)
pam.res2 <- pam(t(matrix_bin_5000bp), 7)
pdf( file = paste(myOutDir_sub2, "13_k-medoids.pdf",  sep="/") )
fviz_cluster(pam.res1)
fviz_cluster(pam.res2)
dev.off() 





# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "complete" )   
res.hc2 <- hclust(res.dist2, method = "complete" )   
res.hc3 <- hclust(res.dist3, method = "complete" )   
res.hc4 <- hclust(res.dist4, method = "complete" )   
pdf( file = paste(myOutDir_sub2, "14A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 



# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D" )   
res.hc2 <- hclust(res.dist2, method = "ward.D" )   
res.hc3 <- hclust(res.dist3, method = "ward.D" )   
res.hc4 <- hclust(res.dist4, method = "ward.D" )   
pdf( file = paste(myOutDir_sub2, "14B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D2" )   
res.hc2 <- hclust(res.dist2, method = "ward.D2" )   
res.hc3 <- hclust(res.dist3, method = "ward.D2" )   
res.hc4 <- hclust(res.dist4, method = "ward.D2" )   
pdf( file = paste(myOutDir_sub2, "14C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




## Fuzzy clustering
res.fanny1 <- fanny( t(matrix_bin_5000bp), 4)
res.fanny2 <- fanny( t(matrix_bin_5000bp), 5)
res.fanny3 <- fanny( t(matrix_bin_5000bp), 6)
res.fanny4 <- fanny( t(matrix_bin_5000bp), 7)

pdf( file = paste(myOutDir_sub2, "15_Fuzzy clustering.pdf",  sep="/") )
fviz_cluster(res.fanny1, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny2, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny3, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny4, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
dev.off() 
























######################################### 
myOutDir_sub2 = paste(outDir, "/6-bin-10000bp",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

myPCA2 <- prcomp( t(matrix_bin_10000bp) )
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




myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = colnames(matrix_1)
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



#####################



dim(matrix_1)
dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)



library("cluster")
library("factoextra")
library("magrittr")


colnames(matrix_1)
colnames(matrix_bin_100bp) <- colnames(matrix_1)
colnames(matrix_bin_500bp) <- colnames(matrix_1)
colnames(matrix_bin_1000bp) <- colnames(matrix_1)
colnames(matrix_bin_5000bp) <- colnames(matrix_1)
colnames(matrix_bin_10000bp) <- colnames(matrix_1)
colnames(matrix_bin_50000bp) <- colnames(matrix_1)


res.dist1 <- get_dist( t(matrix_bin_10000bp) , stand = TRUE, method = "euclidean")
res.dist2 <- get_dist( t(matrix_bin_10000bp) , stand = TRUE, method = "manhattan")
res.dist3 <- get_dist( t(matrix_bin_10000bp) , stand = TRUE, method = "spearman")
res.dist4 <- get_dist( t(matrix_bin_10000bp) , stand = TRUE, method = "pearson")

pdf( file = paste(myOutDir_sub2, "10_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 


pdf( file = paste(myOutDir_sub2, "11_Determining-optimal-number-of-clusters.pdf",  sep="/") )
fviz_nbclust(t(matrix_bin_10000bp), hcut,   method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
fviz_nbclust(t(matrix_bin_10000bp), kmeans, method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
dev.off() 


km.res1 <- kmeans( t(matrix_bin_10000bp), 6, nstart = 25)
km.res2 <- kmeans( t(matrix_bin_10000bp), 7, nstart = 25)

pdf( file = paste(myOutDir_sub2, "12_kmeans.pdf",  sep="/") )
fviz_cluster(km.res1, data =  t(matrix_bin_10000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
fviz_cluster(km.res2, data =  t(matrix_bin_10000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
dev.off() 


## the k-medoids/pam clustering 
pam.res1 <- pam(t(matrix_bin_10000bp), 6)
pam.res2 <- pam(t(matrix_bin_10000bp), 7)
pdf( file = paste(myOutDir_sub2, "13_k-medoids.pdf",  sep="/") )
fviz_cluster(pam.res1)
fviz_cluster(pam.res2)
dev.off() 





# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "complete" )   
res.hc2 <- hclust(res.dist2, method = "complete" )   
res.hc3 <- hclust(res.dist3, method = "complete" )   
res.hc4 <- hclust(res.dist4, method = "complete" )   
pdf( file = paste(myOutDir_sub2, "14A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 



# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D" )   
res.hc2 <- hclust(res.dist2, method = "ward.D" )   
res.hc3 <- hclust(res.dist3, method = "ward.D" )   
res.hc4 <- hclust(res.dist4, method = "ward.D" )   
pdf( file = paste(myOutDir_sub2, "14B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D2" )   
res.hc2 <- hclust(res.dist2, method = "ward.D2" )   
res.hc3 <- hclust(res.dist3, method = "ward.D2" )   
res.hc4 <- hclust(res.dist4, method = "ward.D2" )   
pdf( file = paste(myOutDir_sub2, "14C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




## Fuzzy clustering
res.fanny1 <- fanny( t(matrix_bin_10000bp), 4)
res.fanny2 <- fanny( t(matrix_bin_10000bp), 5)
res.fanny3 <- fanny( t(matrix_bin_10000bp), 6)
res.fanny4 <- fanny( t(matrix_bin_10000bp), 7)

pdf( file = paste(myOutDir_sub2, "15_Fuzzy clustering.pdf",  sep="/") )
fviz_cluster(res.fanny1, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny2, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny3, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny4, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
dev.off() 

























######################################### 
myOutDir_sub2 = paste(outDir, "/7-bin-50000bp",  sep="");
if( ! file.exists(myOutDir_sub2) ) { dir.create(myOutDir_sub2, recursive = TRUE) }

myPCA2 <- prcomp( t(matrix_bin_50000bp) )
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




myPCA2_matrix <- myPCA2$x
dim(myPCA2_matrix)

myLabel = colnames(matrix_1)
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


#####################



dim(matrix_1)
dim(matrix_bin_100bp)
dim(matrix_bin_500bp)
dim(matrix_bin_1000bp)
dim(matrix_bin_5000bp)
dim(matrix_bin_10000bp)
dim(matrix_bin_50000bp)



library("cluster")
library("factoextra")
library("magrittr")


colnames(matrix_1)
colnames(matrix_bin_100bp) <- colnames(matrix_1)
colnames(matrix_bin_500bp) <- colnames(matrix_1)
colnames(matrix_bin_1000bp) <- colnames(matrix_1)
colnames(matrix_bin_5000bp) <- colnames(matrix_1)
colnames(matrix_bin_10000bp) <- colnames(matrix_1)
colnames(matrix_bin_50000bp) <- colnames(matrix_1)


res.dist1 <- get_dist( t(matrix_bin_50000bp) , stand = TRUE, method = "euclidean")
res.dist2 <- get_dist( t(matrix_bin_50000bp) , stand = TRUE, method = "manhattan")
res.dist3 <- get_dist( t(matrix_bin_50000bp) , stand = TRUE, method = "spearman")
res.dist4 <- get_dist( t(matrix_bin_50000bp) , stand = TRUE, method = "pearson")

pdf( file = paste(myOutDir_sub2, "10_visualizing-distance-matrix.pdf",  sep="/") )
fviz_dist(res.dist1,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist2,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist3,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
fviz_dist(res.dist4,  gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")) 
dev.off() 

 
pdf( file = paste(myOutDir_sub2, "11_Determining-optimal-number-of-clusters.pdf",  sep="/") )
fviz_nbclust(t(matrix_bin_50000bp), hcut,   method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
fviz_nbclust(t(matrix_bin_50000bp), kmeans, method = "gap_stat", k.max=10)  ## Determining the optimal number of clusters: 
dev.off() 


km.res1 <- kmeans( t(matrix_bin_50000bp), 6, nstart = 25)
km.res2 <- kmeans( t(matrix_bin_50000bp), 7, nstart = 25)

pdf( file = paste(myOutDir_sub2, "12_kmeans.pdf",  sep="/") )
fviz_cluster(km.res1, data =  t(matrix_bin_50000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
fviz_cluster(km.res2, data =  t(matrix_bin_50000bp),
             ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal())
dev.off() 


## the k-medoids/pam clustering 
pam.res1 <- pam(t(matrix_bin_50000bp), 6)
pam.res2 <- pam(t(matrix_bin_50000bp), 7)
pdf( file = paste(myOutDir_sub2, "13_k-medoids.pdf",  sep="/") )
fviz_cluster(pam.res1)
fviz_cluster(pam.res2)
dev.off() 





# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "complete" )   
res.hc2 <- hclust(res.dist2, method = "complete" )   
res.hc3 <- hclust(res.dist3, method = "complete" )   
res.hc4 <- hclust(res.dist4, method = "complete" )   
pdf( file = paste(myOutDir_sub2, "14A_hierarchical-complete.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 


 
# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D" )   
res.hc2 <- hclust(res.dist2, method = "ward.D" )   
res.hc3 <- hclust(res.dist3, method = "ward.D" )   
res.hc4 <- hclust(res.dist4, method = "ward.D" )   
pdf( file = paste(myOutDir_sub2, "14B_hierarchical-ward.D.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




# Compute hierarchical clustering
res.hc1 <- hclust(res.dist1, method = "ward.D2" )   
res.hc2 <- hclust(res.dist2, method = "ward.D2" )   
res.hc3 <- hclust(res.dist3, method = "ward.D2" )   
res.hc4 <- hclust(res.dist4, method = "ward.D2" )   
pdf( file = paste(myOutDir_sub2, "14C_hierarchical-ward.D2.pdf",  sep="/") )
fviz_dend(res.hc1, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc2, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc3, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
fviz_dend(res.hc4, k = 6,  cex = 0.5,  color_labels_by_k = TRUE,  rect = TRUE )
dev.off() 




## Fuzzy clustering
res.fanny1 <- fanny( t(matrix_bin_50000bp), 4)
res.fanny2 <- fanny( t(matrix_bin_50000bp), 5)
res.fanny3 <- fanny( t(matrix_bin_50000bp), 6)
res.fanny4 <- fanny( t(matrix_bin_50000bp), 7)

pdf( file = paste(myOutDir_sub2, "15_Fuzzy clustering.pdf",  sep="/") )
fviz_cluster(res.fanny1, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny2, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny3, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
fviz_cluster(res.fanny4, ellipse.type = "norm", repel = TRUE,
             palette = "jco", ggtheme = theme_minimal(), legend = "right")
dev.off() 





