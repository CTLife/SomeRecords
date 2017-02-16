
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



MyBoxViolinPlot_1 <- function(vector2,   sampleType2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType) ) + scale_x_discrete(limits=sampleRank2)  +
    geom_errorbar( aes(ymin=min, ymax=max),  data=whisk_1(DataFrame_Local),   width=0.2, size=0.5 ) +
    geom_boxplot( width=0.6,   aes(y=yAxis), outlier.colour="gray45",  outlier.shape=1, outlier.size=0.1, size=0.5, fill=colours2 ) +    
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot",                     sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) +   scale_x_discrete(limits=sampleRank2)  +
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-3adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
    geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "grey30", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-noAdjust",         sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp4 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray") + scale_x_discrete(limits=sampleRank2)  +  
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-noAdjust-colour",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 3) + scale_x_discrete(limits=sampleRank2)  +  
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-3Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
}  






 
MyHistogram_1 <- function(vector2, sampleType2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4 ) {
  dataframeB  <- data.frame( xAxis = vector2,  sampleType=sampleType2 ) 
  
  FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis) )  + geom_bar(aes(fill = sampleType), position = "fill") +  
    xlab(xLab2) + ylab("Relative frequency") +  ggtitle(title2)  +  scale_x_discrete(limits=c("Lowest", "Low", "High", "Highest"))  +
    MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_relative",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType) )  + geom_bar(aes(fill = sampleType) ) +  
    xlab(xLab2) + ylab("Absolute frequency") +  ggtitle(title2)  +  scale_x_discrete(limits=c("Lowest", "Low", "High", "Highest"))  +
    MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_absolute",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
}









WT_rawMatrix1  <- read.table("EEDP5/4-8-1B-WT-LogLinearModel_NTR.txt",   header=FALSE,   sep="",   quote = "",   comment.char = "")  
WT_rawMatrix2  <- read.table("EEDP5_SUZ12/4-8-1B-WT-LogLinearModel_NTR.txt",     header=FALSE,   sep="",   quote = "",   comment.char = "")  
WT_rawMatrix3  <- read.table("SUZ12/4-8-1B-WT-LogLinearModel_NTR.txt",    header=FALSE,   sep="",   quote = "",   comment.char = "")  

Hete_rawMatrix4  <- read.table("EEDP5/4-8-2A-2-EEDheto-LogLinearModel_NTR.txt",   header=FALSE,   sep="",   quote = "",   comment.char = "")  
Hete_rawMatrix5  <- read.table("EEDP5_SUZ12/4-8-2A-2-EEDheto-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  
Hete_rawMatrix6  <- read.table("SUZ12/4-8-2A-2-EEDheto-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  

CKO_rawMatrix7  <- read.table("EEDP5/4-8-2B-2-EEDko-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  
CKO_rawMatrix8  <- read.table("EEDP5_SUZ12/4-8-2B-2-EEDko-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  
CKO_rawMatrix9  <- read.table("SUZ12/4-8-2B-2-EEDko-LogLinearModel_NTR.txt",  header=FALSE,   sep="",   quote = "",   comment.char = "")  






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




####################################################################################################
WT_NTR_merge1   <-  c( WT_NTR1, WT_NTR2,  WT_NTR3  ) 
WT_Label_merge1 <-  c( rep("spe1", length(WT_NTR1)),   rep("common", length(WT_NTR2)),  rep("spe2", length(WT_NTR3)) )  
length(WT_NTR_merge1)
length(WT_Label_merge1)
WT_NTR_merge1[1:100]
WT_Label_merge1[1:100]

index1 <- order(WT_NTR_merge1)        
length(index1)
WT_NTR_merge1[ index1[1:100] ]

WT_NTR_merge2 <- WT_NTR_merge1[index1]
WT_Label_merge2 <- WT_Label_merge1[index1]
length(WT_NTR_merge2)
length(WT_Label_merge2)
WT_NTR_merge2[1:100]
WT_Label_merge2[1:100]

WT_oneClass <- length(WT_NTR_merge2)/4
WT_oneClass2 <- length(WT_NTR_merge2) - 3*WT_oneClass
WT_xAxis <-  c( rep("Lowest", WT_oneClass),   rep("Low", WT_oneClass),  rep("High", WT_oneClass) ,  rep("Highest", WT_oneClass2+2 ) ) 
length(WT_xAxis)




MyHistogram_1(vector2=WT_xAxis, sampleType2=WT_Label_merge2,  path2="z-figures",   
              fileName2="5-WT-histogram",  title2="4 classes of NTR",  xLab2="NTR", height2=4,  width2=4 )
  
####################################################################################################


 










####################################################################################################
Hete_NTR_merge1   <-  c( Hete_NTR4, Hete_NTR5,  Hete_NTR6  ) 
Hete_Label_merge1 <-  c( rep("spe1", length(Hete_NTR4)),   rep("common", length(Hete_NTR5)),  rep("spe2", length(Hete_NTR6)) )  
length(Hete_NTR_merge1)
length(Hete_Label_merge1)
Hete_NTR_merge1[1:100]
Hete_Label_merge1[1:100]

index1 <- order(Hete_NTR_merge1)        
length(index1)
Hete_NTR_merge1[ index1[1:100] ]

Hete_NTR_merge2 <- Hete_NTR_merge1[index1]
Hete_Label_merge2 <- Hete_Label_merge1[index1]
length(Hete_NTR_merge2)
length(Hete_Label_merge2)
Hete_NTR_merge2[1:100]
Hete_Label_merge2[1:100]

Hete_oneClass <- length(Hete_NTR_merge2)/4
Hete_oneClass2 <- length(Hete_NTR_merge2) - 3*Hete_oneClass
Hete_xAxis <-  c( rep("Lowest", Hete_oneClass),   rep("Low", Hete_oneClass),  rep("High", Hete_oneClass) ,  rep("Highest", Hete_oneClass2+2 ) ) 
length(Hete_xAxis)




MyHistogram_1(vector2=Hete_xAxis, sampleType2=Hete_Label_merge2,  path2="z-figures",   
              fileName2="5-Hete-histogram",  title2="4 classes of NTR",  xLab2="NTR", height2=4,  width2=4 )

####################################################################################################











####################################################################################################
CKO_NTR_merge1   <-  c( CKO_NTR7, CKO_NTR8,  CKO_NTR9  ) 
CKO_Label_merge1 <-  c( rep("spe1", length(CKO_NTR7)),   rep("common", length(CKO_NTR8)),  rep("spe2", length(CKO_NTR9)) )  
length(CKO_NTR_merge1)
length(CKO_Label_merge1)
CKO_NTR_merge1[1:100]
CKO_Label_merge1[1:100]

index1 <- order(CKO_NTR_merge1)        
length(index1)
CKO_NTR_merge1[ index1[1:100] ]

CKO_NTR_merge2 <- CKO_NTR_merge1[index1]
CKO_Label_merge2 <- CKO_Label_merge1[index1]
length(CKO_NTR_merge2)
length(CKO_Label_merge2)
CKO_NTR_merge2[1:100]
CKO_Label_merge2[1:100]

CKO_oneClass <- length(CKO_NTR_merge2)/4
CKO_oneClass2 <- length(CKO_NTR_merge2) - 3*CKO_oneClass
CKO_xAxis <-  c( rep("Lowest", CKO_oneClass),   rep("Low", CKO_oneClass),  rep("High", CKO_oneClass) ,  rep("Highest", CKO_oneClass2+2) ) 
length(CKO_xAxis)




MyHistogram_1(vector2=CKO_xAxis, sampleType2=CKO_Label_merge2,  path2="z-figures",   
              fileName2="5-CKO-histogram",  title2="4 classes of NTR",  xLab2="NTR", height2=4,  width2=4 )

####################################################################################################





















length( WT_NTR1 )
length( WT_NTR2 )
length( WT_NTR3 )
length( Hete_NTR4 )
length( Hete_NTR5 )
length( Hete_NTR6 )
length( CKO_NTR7 )
length( CKO_NTR8 )
length( CKO_NTR9 )






####################################################################################################
WT_NTR_merge1   <-  c( WT_NTR1[1:7000], WT_NTR2,  WT_NTR3[1:8000]  ) 
WT_Label_merge1 <-  c( rep("spe1", length(WT_NTR1[1:7000])),   rep("common", length(WT_NTR2)),  rep("spe2", length(WT_NTR3[1:8000])) )  
length(WT_NTR_merge1)
length(WT_Label_merge1)
WT_NTR_merge1[1:100]
WT_Label_merge1[1:100]

index1 <- order(WT_NTR_merge1)        
length(index1)
WT_NTR_merge1[ index1[1:100] ]

WT_NTR_merge2 <- WT_NTR_merge1[index1]
WT_Label_merge2 <- WT_Label_merge1[index1]
length(WT_NTR_merge2)
length(WT_Label_merge2)
WT_NTR_merge2[1:100]
WT_Label_merge2[1:100]

WT_oneClass <- length(WT_NTR_merge2)/4
WT_oneClass2 <- length(WT_NTR_merge2) - 3*WT_oneClass
WT_xAxis <-  c( rep("Lowest", WT_oneClass),   rep("Low", WT_oneClass),  rep("High", WT_oneClass) ,  rep("Highest", WT_oneClass2+1 ) ) 
length(WT_xAxis)




MyHistogram_1(vector2=WT_xAxis, sampleType2=WT_Label_merge2,  path2="z-figures",   
              fileName2="6-WT-histogram",  title2="4 classes of NTR",  xLab2="NTR", height2=4,  width2=4 )

####################################################################################################













####################################################################################################
Hete_NTR_merge1   <-  c( Hete_NTR4[1:7000], Hete_NTR5,  Hete_NTR6[1:8000]  ) 
Hete_Label_merge1 <-  c( rep("spe1", length(Hete_NTR4[1:7000])),   rep("common", length(Hete_NTR5)),  rep("spe2", length(Hete_NTR6[1:8000])) )  
length(Hete_NTR_merge1)
length(Hete_Label_merge1)
Hete_NTR_merge1[1:100]
Hete_Label_merge1[1:100]

index1 <- order(Hete_NTR_merge1)        
length(index1)
Hete_NTR_merge1[ index1[1:100] ]

Hete_NTR_merge2 <- Hete_NTR_merge1[index1]
Hete_Label_merge2 <- Hete_Label_merge1[index1]
length(Hete_NTR_merge2)
length(Hete_Label_merge2)
Hete_NTR_merge2[1:100]
Hete_Label_merge2[1:100]

Hete_oneClass <- length(Hete_NTR_merge2)/4
Hete_oneClass2 <- length(Hete_NTR_merge2) - 3*Hete_oneClass
Hete_xAxis <-  c( rep("Lowest", Hete_oneClass),   rep("Low", Hete_oneClass),  rep("High", Hete_oneClass) ,  rep("Highest", Hete_oneClass2+1 ) ) 
length(Hete_xAxis)




MyHistogram_1(vector2=Hete_xAxis, sampleType2=Hete_Label_merge2,  path2="z-figures",   
              fileName2="6-Hete-histogram",  title2="4 classes of NTR",  xLab2="NTR", height2=4,  width2=4 )

####################################################################################################











####################################################################################################
CKO_NTR_merge1   <-  c( CKO_NTR7[1:7000], CKO_NTR8,  CKO_NTR9[1:8000]  ) 
CKO_Label_merge1 <-  c( rep("spe1", length(CKO_NTR7[1:7000])),   rep("common", length(CKO_NTR8)),  rep("spe2", length(CKO_NTR9[1:8000])) )  
length(CKO_NTR_merge1)
length(CKO_Label_merge1)
CKO_NTR_merge1[1:100]
CKO_Label_merge1[1:100]

index1 <- order(CKO_NTR_merge1)        
length(index1)
CKO_NTR_merge1[ index1[1:100] ]

CKO_NTR_merge2 <- CKO_NTR_merge1[index1]
CKO_Label_merge2 <- CKO_Label_merge1[index1]
length(CKO_NTR_merge2)
length(CKO_Label_merge2)
CKO_NTR_merge2[1:100]
CKO_Label_merge2[1:100]

CKO_oneClass <- length(CKO_NTR_merge2)/4
CKO_oneClass2 <- length(CKO_NTR_merge2) - 3*CKO_oneClass
CKO_xAxis <-  c( rep("Lowest", CKO_oneClass),   rep("Low", CKO_oneClass),  rep("High", CKO_oneClass) ,  rep("Highest", CKO_oneClass2+1) ) 
length(CKO_xAxis)




MyHistogram_1(vector2=CKO_xAxis, sampleType2=CKO_Label_merge2,  path2="z-figures",   
              fileName2="6-CKO-histogram",  title2="4 classes of NTR",  xLab2="NTR", height2=4,  width2=4 )

####################################################################################################



















