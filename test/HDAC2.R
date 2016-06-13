

AllResults_g <- "HDAC2_figures"
if( ! file.exists(AllResults_g) ) { dir.create(AllResults_g) }

library("reshape2")
library("ggplot2") 
library("grid")
library("Cairo")
library("RColorBrewer")
library("gplots")  
library("stats")
library("KernSmooth")
library("psych")
library("minerva")
library("matrixStats")
library("coin")
library("svglite")
library("extrafont")
font_import()
fonts()
fonttable()
loadfonts()
loadfonts(device="postscript")
names(postscriptFonts())


MyTheme_1 <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    #hjust=1,    vjust=1,      angle=30
        theme(  
                line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                              ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
                rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                       ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
                text  = element_text(family="serif",  face=NULL,  colour="black",  size=textSize1, hjust=NULL, vjust=NULL,   angle=NULL, lineheight=NULL),     ## all text elements.           "serif" for a serif font. 所有文本相关属性.
                title = element_text(family="serif",  face=NULL,  colour="black",  size=textSize1, hjust=NULL, vjust=NULL,   angle=NULL, lineheight=NULL),     ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    
                axis.title        = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=NULL,   angle=NULL,   lineheight=NULL  ),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
                axis.title.x      = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=-0.5,   angle=NULL,   lineheight=NULL  ),       ## x axis label (element_text; inherits from axis.title)
                axis.title.y      = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=1.5,    angle=NULL,   lineheight=NULL  ),       ## y axis label (element_text; inherits from axis.title)
                axis.text         = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=NULL,   angle=NULL,   lineheight=NULL  ),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
                axis.text.x       = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1, lineheight=angle1),       ## x axis tick labels (element_text; inherits from axis.text)
                axis.text.y       = element_text(family="serif", face=NULL, colour="black", size=textSize1,    hjust=NULL,   vjust=NULL,   angle=NULL,   lineheight=NULL  ),       ## y axis tick labels (element_text; inherits from axis.text)
                axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),                                                                              ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
                axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),                                                                              ## x axis tick marks (element_line; inherits from axis.ticks)
                axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),                                                                              ## y axis tick marks (element_line; inherits from axis.ticks)
                axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                                                                                          ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
                #axis.text.margin = grid::unit(1.0,   "mm",   data=NULL),  	                                                                                                   ## space between tick mark and tick label (unit),  ‘"mm"’ Millimetres.  10 mm = 1 cm. 刻度线和刻度标签之间的间距.                                                                           
                axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                                   ## lines along axes (element_line; inherits from line). 坐标轴线
                axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                                   ## line along x axis (element_line; inherits from axis.line)
                axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	                                                                   ## line along y axis (element_line; inherits from axis.line)
    
                legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
                legend.margin        = grid::unit(1, "mm", data=NULL), 	                                                      ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
                legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
                legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                              ## size of legend keys   (unit; inherits from legend.key.size)
                legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                              ## key background height (unit; inherits from legend.key.size)
                legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                     ## key background width  (unit; inherits from legend.key.size)
                legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
                legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
                legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
                legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
                legend.position      = "right", 	            ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
                legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
                legend.justification = "center",      	            ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
                legend.box           = NULL, 	                    ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
                legend.box.just      = NULL, 	                    ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
                panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
                panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
                panel.margin       = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
                panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
                panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	## major grid lines (element_line; inherits from panel.grid)  主网格线
                panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
                panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	## vertical major grid lines (element_line; inherits from panel.grid.major)
                panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
                panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
                panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
                plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
                plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=1.5,   angle=NULL, lineheight=NULL),    ## plot title (text appearance) (element_text; inherits from title)  图形标题
                plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
                strip.background = element_rect(colour=NULL, size=NULL, linetype=NULL, fill=NULL ), 	                                                        ## background of facet labels (element_rect; inherits from rect)  分面标签背景
                strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	## facet labels (element_text; inherits from text)
                strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	## facet labels along horizontal direction (element_text; inherits from strip.text)
                strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	## facet labels along vertical direction (element_text; inherits from strip.text) 
        ) 
} 







MySaveGgplot2_1 <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
        SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
        PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
        PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
        EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
        EPS2 <- paste(path1,  "/",  "EPS2", sep = "",  collapse = NULL)
        if( ! file.exists(SVG1) ) { dir.create(SVG1) }
        if( ! file.exists(PNG1) ) { dir.create(PNG1) }
        if( ! file.exists(PDF1) ) { dir.create(PDF1) }
        if( ! file.exists(EPS1) ) { dir.create(EPS1) }
        if( ! file.exists(EPS2) ) { dir.create(EPS2) }
        postscript( file=paste(EPS2,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),   height =height1, width = width1,  family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)      
        print( ggplot2Figure1 )
        dev.off()  
        ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000 )
        ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1000,   device=cairo_ps)         
}
                            
  

   


MyAverageLines_1 <- function(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
        binLen <- 20    ## Bin size is 20 bp
        binNum <- 500  ## 500 * 20 = 10000 bp
        Position_Local   <-  seq(from = -5,  by=0.02,  length.out=binNum)   ## unit is kb
        DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2 ) 
  
        FigureTemp1 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                        scale_colour_manual( values=colours2, breaks=sampleRank2  ) + geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                        scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5", center2,  "2.5",  "5") ) +  
                        MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,  height1=height2, width1=width2)
        
        FigureTemp2 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  + xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2 ) +
                        scale_colour_manual( values=colours2, breaks=sampleRank2  ) +  geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                        scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5",  center2,  "2.5",  "5") ) +  
                        MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  



MyAverageLines_2 <- function(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
  binLen <- 20    ## Bin size is 20 bp
  binNum <- 100  ## 500 * 20 = 10000 bp
  Position_Local   <-  seq(from = -1,  by=0.02,  length.out=binNum)   ## unit is kb
  DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2 ) 
  
  FigureTemp1 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  +  xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
    scale_colour_manual( values=colours2, breaks=sampleRank2  ) + geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
    scale_x_continuous( breaks=c( -1,   -0.5,  0,  0.5,  1  ), labels=c( "-1",   "-0.5",  center2,  "0.5",  "1" ) ) +  
    MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot(DataFrame_Local,   aes(x=xAxis, y=yAxis,  color=SampleType) )  + xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + ylim(Ymin2, Ymax2 ) +
    scale_colour_manual( values=colours2, breaks=sampleRank2  ) +  geom_line(size=0.5) +  geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
    scale_x_continuous( breaks=c( -1,   -0.5,  0,  0.5,  1  ), labels=c( "-1",   "-0.5",  center2,  "0.5",  "1"  ) ) +  
    MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  





normalize_1 <- function(vector1) {
  max1 = max(vector1)  
  min1 = min(vector1)  
  lower1 = -1    
  upper1 = 1
  vector2 = lower1 + (upper1-lower1)*(vector1-min1)/(max1-min1)
  return(vector2)
}



normalize_2 <- function(matrix1) {
  matrix2 <- matrix1 
  numc <- nrow(matrix1)
  for(i in c(1:numc)) {matrix2[i,] <- normalize_1(matrix1[i,]) }
  return(matrix2)
}





#######################################################################################################################################################################
down_HetoHDAC2   <- read.table("102downGenes-EED-P5peaks/heto-lib1644-HDAC2_center_heatmap/102downGenes-EED-P5peaks.bed.1-BED_1_heto-lib1644-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  
down_HomoHDAC2   <- read.table("102downGenes-EED-P5peaks/homo-lib1642-HDAC2_center_heatmap/102downGenes-EED-P5peaks.bed.1-BED_1_homo-lib1642-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  

up_HetoHDAC2   <- read.table("303upGenes-EED-P5peaks/heto-lib1644-HDAC2_center_heatmap/303upGenes-EED-P5peaks.bed.1-BED_1_heto-lib1644-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  
up_HomoHDAC2   <- read.table("303upGenes-EED-P5peaks/homo-lib1642-HDAC2_center_heatmap/303upGenes-EED-P5peaks.bed.1-BED_1_homo-lib1642-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  

retain_HetoHDAC2   <- read.table("9631retainGenes-EED-P5peaks/heto-lib1644-HDAC2_center_heatmap/9631retainGenes-EED-P5peaks.bed.1-BED_1_heto-lib1644-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  
retain_HomoHDAC2   <- read.table("9631retainGenes-EED-P5peaks/homo-lib1642-HDAC2_center_heatmap/9631retainGenes-EED-P5peaks.bed.1-BED_1_homo-lib1642-HDAC2-adultCMs_Rep1.bgsub.Fnor.smooth.wig.heatmap.xls",                 header=TRUE,   sep="",   quote = "",   comment.char = "")  


dim(down_HetoHDAC2)
dim(down_HomoHDAC2)
dim(up_HetoHDAC2)
dim(up_HomoHDAC2)
dim(retain_HetoHDAC2)
dim(retain_HomoHDAC2)


down_HetoHDAC2   <- as.matrix(down_HetoHDAC2[, -c(1:4)])
down_HomoHDAC2   <- as.matrix(down_HomoHDAC2[, -c(1:4)])
up_HetoHDAC2   <- as.matrix(up_HetoHDAC2[, -c(1:4)])
up_HomoHDAC2   <- as.matrix(up_HomoHDAC2[, -c(1:4)])
retain_HetoHDAC2   <- as.matrix(retain_HetoHDAC2[, -c(1:4)])
retain_HomoHDAC2   <- as.matrix(retain_HomoHDAC2[, -c(1:4)])


dim(down_HetoHDAC2)
dim(down_HomoHDAC2)
dim(up_HetoHDAC2)
dim(up_HomoHDAC2)
dim(retain_HetoHDAC2)
dim(retain_HomoHDAC2)


sink(  file=paste(AllResults_g, "/", "summary.txt",  sep="")  )
summary( as.vector(down_HetoHDAC2) )
summary( as.vector(down_HomoHDAC2) )
summary( as.vector(up_HetoHDAC2) )
summary( as.vector(up_HomoHDAC2) )
summary( as.vector(retain_HetoHDAC2) )
summary( as.vector(retain_HomoHDAC2) )
sink()

down_HetoHDAC2A   <- as.matrix(down_HetoHDAC2)
down_HomoHDAC2A   <- as.matrix(down_HomoHDAC2)
up_HetoHDAC2A   <- as.matrix(up_HetoHDAC2)
up_HomoHDAC2A   <- as.matrix(up_HomoHDAC2)
retain_HetoHDAC2A   <- as.matrix(retain_HetoHDAC2)
retain_HomoHDAC2A   <- as.matrix(retain_HomoHDAC2)

dim(down_HetoHDAC2A)
dim(down_HomoHDAC2A)
dim(up_HetoHDAC2A)
dim(up_HomoHDAC2A)
dim(retain_HetoHDAC2A)
dim(retain_HomoHDAC2A)

numOfColumns2 <- ncol(down_HetoHDAC2A)
numOfColumns2




MyAverageLines_3 <- function(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
        binLen <- 20    ## Bin size is 20 bp
        binNum <- 500  ## 500 * 20 = 10000 bp
        Position_Local   <-  seq(from = -5,  by=0.02,  length.out=binNum)   ## unit is kb
        DataFrame_Local  <-  data.frame( xAxis = c( rep(Position_Local, numSample2) ),      yAxis = vector2,    SampleType = sampleType2 ) 
  
        FigureTemp1 <- ggplot(DataFrame_Local )  +  
                        geom_line(  aes(x=xAxis,  y=yAxis,  color=SampleType, linetype=SampleType, size=SampleType)  )+ 
                        xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                        scale_colour_manual( values=colours2   ) +
                        scale_linetype_manual(  values=c(1, 3,  1, 3, 1, 3 )  ) + 
                        scale_size_manual(values = c(0.2, 0.4, 0.2, 0.4, 0.2, 0.4)  )+
                        geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   
                        scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5", center2,  "2.5",  "5") ) +  
                        MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=fileName2,  height1=height2, width1=width2)
        
        FigureTemp2 <- ggplot(DataFrame_Local  )  +  
                        geom_line(  aes(x=xAxis,  y=yAxis,  color=SampleType, linetype=SampleType, size=SampleType)  )+ 
                        xlab(xLab2) +  ylab(yLab2) +  ggtitle(title2) + 
                        scale_colour_manual( values=colours2   ) +
                        scale_linetype_manual(  values=c(1, 3,  1, 3, 1, 3 )  ) + 
                        scale_size_manual(values = c(0.2, 0.4, 0.2, 0.4, 0.2, 0.4)  )  + 
                        geom_vline(xintercept=0, lty=2, col="gray45", size=0.5) +   ylim(Ymin2, Ymax2 ) +
                        scale_x_continuous( breaks=c(-5,  -2.5,  0,  2.5,   5  ), labels=c("-5",  "-2.5",  center2,  "2.5",  "5") ) +  
                        MyTheme_1(textSize1=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL ) +  guides( colour = guide_legend(override.aes = list(size=1.5)) )  
        MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "-limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  




MyAverageLines_3(vector2=c( colMeans(down_HetoHDAC2A ),  colMeans(down_HomoHDAC2A), colMeans(up_HetoHDAC2A ),  colMeans(up_HomoHDAC2A),
                            colMeans(retain_HetoHDAC2A ),  colMeans(retain_HomoHDAC2A) ), 
                 numSample2=6,   
                 sampleType2=c( rep("down_HDAC2_Heto", numOfColumns2),      rep("down_HDAC2_Homo", numOfColumns2), rep("up_HDAC2_Heto", numOfColumns2),    rep("up_HDAC2_Homo", numOfColumns2), 
                                rep("unchange_HDAC2_Heto", numOfColumns2),    rep("unchange_HDAC2_Homo", numOfColumns2)  ), 
                 sampleRank2=c( "down_HDAC2_Heto",  "down_HDAC2_Homo",    "up_HDAC2_Heto",  "up_HDAC2_Homo",   "unchange_HDAC2_Heto",  "unchange_HDAC2_Homo"   ),     
                 colours2=c( "down_HDAC2_Heto"="blue",  "down_HDAC2_Homo"="blue",    "up_HDAC2_Heto"="red",  "up_HDAC2_Homo"="red",   
                             "unchange_HDAC2_Heto"="black",  "unchange_HDAC2_Homo"="black"  ), 
                 path2=AllResults_g,     fileName2=paste(AllResults_g, "_averageCurve", sep=""),  
                 title2="EED peaks is associated with DEGs",     xLab2="Relative distance (kb)",    yLab2="HDAC2 reads density",   
                 Ymin2=0.3,   Ymax2=3.5,    height2=3.5,   width2=7,   center2="0" )




MyAverageLines_3(vector2=c( colMeans(down_HetoHDAC2A ),  colMeans(down_HomoHDAC2A) ), 
                 numSample2=2,   
                 sampleType2=c( rep("down_HDAC2_Heto", numOfColumns2),      rep("down_HDAC2_Homo", numOfColumns2)  ), 
                 sampleRank2=c( "down_HDAC2_Heto",  "down_HDAC2_Homo"   ),     
                 colours2=c( "down_HDAC2_Heto"="blue",  "down_HDAC2_Homo"="blue"  ), 
                 path2=AllResults_g,     fileName2=paste(AllResults_g, "_down", sep=""),  
                 title2="EED peaks is associated with down-regulated genes",     xLab2="Relative distance (kb)",    yLab2="HDAC2 reads density",   
                 Ymin2=0.3,   Ymax2=3.5,    height2=3.5,   width2=7,   center2="0" )



MyAverageLines_3(vector2=c(  colMeans(up_HetoHDAC2A ),  colMeans(up_HomoHDAC2A) ), 
                 numSample2=2,   
                 sampleType2=c(  rep("up_HDAC2_Heto", numOfColumns2),    rep("up_HDAC2_Homo", numOfColumns2)  ), 
                 sampleRank2=c(    "up_HDAC2_Heto",  "up_HDAC2_Homo"   ),     
                 colours2=c( "up_HDAC2_Heto"="red",  "up_HDAC2_Homo"="red"  ), 
                 path2=AllResults_g,     fileName2=paste(AllResults_g, "_up", sep=""),  
                 title2="EED peaks is associated with up-regulated genes",     xLab2="Relative distance (kb)",    yLab2="HDAC2 reads density",   
                 Ymin2=0.3,   Ymax2=3.5,    height2=3.5,   width2=7,   center2="0" )



MyAverageLines_3(vector2=c(  colMeans(retain_HetoHDAC2A ),  colMeans(retain_HomoHDAC2A) ), 
                 numSample2=2,   
                 sampleType2=c(  rep("unchange_HDAC2_Heto", numOfColumns2),    rep("unchange_HDAC2_Homo", numOfColumns2)  ), 
                 sampleRank2=c(   "unchange_HDAC2_Heto",  "unchange_HDAC2_Homo"   ),     
                 colours2=c(  "unchange_HDAC2_Heto"="black",  "unchange_HDAC2_Homo"="black"  ), 
                 path2=AllResults_g,     fileName2=paste(AllResults_g, "_unchange", sep=""),  
                 title2="EED peaks is associated with unchanged genes",     xLab2="Relative distance (kb)",    yLab2="HDAC2 reads density",   
                 Ymin2=0.3,   Ymax2=3.5,    height2=3.5,   width2=7,   center2="0" )




