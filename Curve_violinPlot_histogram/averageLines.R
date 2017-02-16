
library("ggplot2") 
library("reshape2") 
library("psych")
library("minerva")

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



## for center±5kb 
MyAverageLines_1 <- function(vector2,  numSample2,  sampleType2, sampleRank2,   colours2,  path2,  fileName2,  title2,  xLab2,  yLab2,   Ymin2=0,   Ymax2=10,   height2=3, width2=4, center2="0") {  
  binLen <- 20    ## Bin size is 20 bp
  binNum <- 500   ## 500 * 20 = 10000 bp
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
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_limitY",  sep="",  collapse=NULL),  height1=height2, width1=width2)
}  




WT_rawMatrix1  <- read.table("spe_H3K27ac/4-1-1B-LogLinearModel_NTR.txt",   header=FALSE,   sep="",   quote = "",   comment.char = "")  
WT_rawMatrix2  <- read.table("common/4-1-1B-LogLinearModel_NTR.txt",     header=FALSE,   sep="",   quote = "",   comment.char = "")  
WT_rawMatrix3  <- read.table("spe_HDAC2/4-1-1B-LogLinearModel_NTR.txt",    header=FALSE,   sep="",   quote = "",   comment.char = "")  

Hete_rawMatrix4  <- read.table("spe_H3K27ac/4-1-2A-LogLinearModel_NTR.txt",   header=FALSE,   sep="",   quote = "",   comment.char = "")  
Hete_rawMatrix5  <- read.table("common/4-1-2A-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  
Hete_rawMatrix6  <- read.table("spe_HDAC2/4-1-2A-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  

CKO_rawMatrix7  <- read.table("spe_H3K27ac/4-1-2B-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  
CKO_rawMatrix8  <- read.table("common/4-1-2B-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  
CKO_rawMatrix9  <- read.table("spe_HDAC2/4-1-2B-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  



dim(WT_rawMatrix1)
dim(WT_rawMatrix2)
dim(WT_rawMatrix3)
dim(Hete_rawMatrix4)
dim(Hete_rawMatrix5)
dim(Hete_rawMatrix6)
dim(CKO_rawMatrix7)
dim(CKO_rawMatrix8)
dim(CKO_rawMatrix9)




WT_NTR1  <- WT_rawMatrix1[,1]
WT_NTR2  <- WT_rawMatrix2[,1]
WT_NTR3  <- WT_rawMatrix3[,1]
Hete_NTR4  <- Hete_rawMatrix4[,1]
Hete_NTR5  <- Hete_rawMatrix5[,1]
Hete_NTR6  <- Hete_rawMatrix6[,1]
CKO_NTR7  <- CKO_rawMatrix7[,1]
CKO_NTR8  <- CKO_rawMatrix8[,1]
CKO_NTR9  <- CKO_rawMatrix9[,1]



length( WT_NTR1 )
length( WT_NTR2 )
length( WT_NTR3 )
length( Hete_NTR4 )
length( Hete_NTR5 )
length( Hete_NTR6 )
length( CKO_NTR7 )
length( CKO_NTR8 )
length( CKO_NTR9 )



MyAverageLines_1(vector2= c( WT_NTR1, WT_NTR2,  WT_NTR3  ),   
                  numSample2=3,   
                  sampleType2=c( rep("WT1", 500),   rep("WT2", 500),  rep("WT3", 500)  ), 
                  
                  sampleRank2=c(  "WT1", "WT2",  "WT3"  ), 
                  
                  colours2=c(  "WT1"="red", "WT2"="blue",  "WT3"="purple"   ), 
                  path2="z-figures",   fileName2= paste("3-WT-averageLines",   "8classes",  sep = "_") ,  
                  title2="NTR",   xLab2="Peaks",    yLab2="NTR",   
                  Ymin2=0,   Ymax2=0.3,    height2=3.25,   width2=4.8, center2="0" )    ## height=4cm, width=5cm






MyAverageLines_1(vector2= c( Hete_NTR4, Hete_NTR5,  Hete_NTR6  ),   
                 numSample2=3,   
                 sampleType2=c( rep("WT1", 500),   rep("WT2", 500),  rep("WT3", 500)  ), 
                 
                 sampleRank2=c(  "WT1", "WT2",  "WT3"  ), 
                 
                 colours2=c(  "WT1"="red", "WT2"="blue",  "WT3"="purple"   ), 
                 path2="z-figures",   fileName2= paste("3-Hete-averageLines",   "8classes",  sep = "_") ,  
                 title2="NTR",   xLab2="Peaks",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.25,   width2=4.8, center2="0" )    ## height=4cm, width=5cm





MyAverageLines_1(vector2= c( CKO_NTR7, CKO_NTR8,  CKO_NTR9  ),   
                 numSample2=3,   
                 sampleType2=c( rep("WT1", 500),   rep("WT2", 500),  rep("WT3", 500)  ), 
                 
                 sampleRank2=c(  "WT1", "WT2",  "WT3"  ), 
                 
                 colours2=c(  "WT1"="red", "WT2"="blue",  "WT3"="purple"   ), 
                 path2="z-figures",   fileName2= paste("3-CKO-averageLines",   "8classes",  sep = "_") ,  
                 title2="NTR",   xLab2="Peaks",    yLab2="NTR",   
                 Ymin2=0,   Ymax2=0.3,    height2=3.25,   width2=4.8, center2="0" )    ## height=4cm, width=5cm















