
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

  FigureTemp5 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 2) + scale_x_discrete(limits=sampleRank2)  +  
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-2Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  
}  



MyGroupedBoxplot_1 <- function(vector2,   sampleType2,   fill2, sampleRank2,   colours2,   path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2, fillcolor = fill2    ) 
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=fillcolor, colour=fillcolor) ) + scale_x_discrete(limits=sampleRank2)  +
    geom_boxplot( position=position_dodge(0.9), outlier.size=0 , alpha=0  ) +    
    geom_jitter( aes(colour=type1 ), size=0.1,  alpha=0.25, position = position_jitterdodge(jitter.width=0.7, dodge.width=0.9 )  ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_groupBox",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
}  



MyGroupedBoxplot_1 <- function(vector2,   sampleType2,   fill2, sampleRank2,  path2,   fileName2,    title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2, fillcolor = fill2    ) 
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=fillcolor, colour=fillcolor) ) + scale_x_discrete(limits=sampleRank2)  +
    geom_boxplot( position=position_dodge(0.9), outlier.size=0 , alpha=0  ) +    
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_groupBox",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp2 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=fillcolor, colour=fillcolor) ) + scale_x_discrete(limits=sampleRank2)  +
    geom_boxplot( position=position_dodge(0.9), outlier.size=0 , alpha=0  ) +    
    geom_jitter( aes(colour=fillcolor ), size=0.1,  alpha=0.25, position = position_jitterdodge(jitter.width=0.7, dodge.width=0.9 )  ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_groupBox_scatter",  sep="",  collapse=NULL),  height1=height2, width1=width2)
  
}  




MyHistogram_1 <- function(vector2, sampleType2,  path2,   fileName2,  title2,  xLab2, height2=4,  width2=4 ) {
  dataframeB  <- data.frame( xAxis = vector2,  sampleType=sampleType2 ) 
  
  FigureTemp1 <- ggplot(data=dataframeB, aes(x=xAxis) )  + geom_bar(aes(fill = sampleType), position = "fill"  ) +  
    scale_fill_manual(values = c("spe1"="purple2", "common"="blue2",  "spe2"="red2") ) + 
    xlab(xLab2) + ylab("Relative frequency") +  ggtitle(title2)  +  scale_x_discrete(limits=c("Lowest", "Low", "High", "Highest"))  +
    MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_relative",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
  
  FigureTemp2 <- ggplot(data=dataframeB, aes(x=xAxis, fill=sampleType) )  + geom_bar(aes(fill = sampleType) ) +  
    scale_fill_manual(values = c("spe1"="purple2", "common"="blue2",  "spe2"="red2") ) + 
    xlab(xLab2) + ylab("Absolute frequency") +  ggtitle(title2)  +  scale_x_discrete(limits=c("Lowest", "Low", "High", "Highest"))  +
    MyTheme_1( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_absolute",    sep="",  collapse=NULL),  height1=height2,  width1=width2)
}





###########################################################################################################################################################################
Housekeeping_1    <- read.table("Housekeeping_NTR_TPM.2303genes",  header=TRUE,   sep="",   quote = "",   comment.char = "")  
Specific_2        <- read.table("Specific_NTR_TPM.1450genes",      header=TRUE,   sep="",   quote = "",   comment.char = "")  
            
dim(Housekeeping_1)
dim(Specific_2)

Housekeeping_1[1:10,]
Specific_2[1:10,]

AllResults_g <- "Z-FinalFigures"
if( ! file.exists(AllResults_g) ) { dir.create(path=AllResults_g,   recursive = TRUE) }



MyBoxViolinPlot_1(vector2= log2( c( Housekeeping_1[,3],   Specific_2[,3]  ) + 1 ),   
                  sampleType2=c( rep("Housekeeping", length(Housekeeping_1[,3])),   rep("HeartSpecific", length(Specific_2[,3]))   ), 
                  sampleRank2=c(  "Housekeeping",   "HeartSpecific"  ), 
                  colours2=c(   "Housekeeping"="red2",    "HeartSpecific"="blue2"   ), 
                  path2=AllResults_g,   fileName2= paste("1-BoxViolin-all",   "TPM",  sep = "_") ,  
                  title2="Genes",   xLab2="Type",    yLab2="log2(TPM+1)",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=12)    ## width = 1 + 4*0.5, height=5cm 



MyBoxViolinPlot_1(vector2=  c( Housekeeping_1[,4],   Specific_2[,4]  ) ,   
                  sampleType2=c( rep("Housekeeping", length(Housekeeping_1[,4])),   rep("HeartSpecific", length(Specific_2[,4]))   ), 
                  sampleRank2=c(  "Housekeeping",   "HeartSpecific"  ), 
                  colours2=c(   "Housekeeping"="red2",    "HeartSpecific"="blue2"   ), 
                  path2=AllResults_g,   fileName2= paste("2-BoxViolin-all",   "NTR-1kb",  sep = "_") ,  
                  title2="Genes",   xLab2="Type",    yLab2="NTR",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=1.0)    ## width = 1 + 4*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=  c( Housekeeping_1[,5],   Specific_2[,5]  ) ,   
                  sampleType2=c( rep("Housekeeping", length(Housekeeping_1[,5])),   rep("HeartSpecific", length(Specific_2[,5]))   ), 
                  sampleRank2=c(  "Housekeeping",   "HeartSpecific"  ), 
                  colours2=c(   "Housekeeping"="red2",    "HeartSpecific"="blue2"   ), 
                  path2=AllResults_g,   fileName2= paste("3-BoxViolin-all",   "NTR-500bp",  sep = "_") ,  
                  title2="Genes",   xLab2="Type",    yLab2="NTR",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=1.0)    ## width = 1 + 4*0.5, height=5cm 





MyBoxViolinPlot_1(vector2=  c( Housekeeping_1[,6],   Specific_2[,6]  ) ,   
                  sampleType2=c( rep("Housekeeping", length(Housekeeping_1[,6])),   rep("HeartSpecific", length(Specific_2[,6]))   ), 
                  sampleRank2=c(  "Housekeeping",   "HeartSpecific"  ), 
                  colours2=c(   "Housekeeping"="red2",    "HeartSpecific"="blue2"   ), 
                  path2=AllResults_g,   fileName2= paste("4-BoxViolin-all",   "NTR-max",  sep = "_") ,  
                  title2="Genes",   xLab2="Type",    yLab2="NTR",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=1.0)    ## width = 1 + 4*0.5, height=5cm 




MyBoxViolinPlot_1(vector2=  c( Housekeeping_1[,7],   Specific_2[,7]  ) ,   
                  sampleType2=c( rep("Housekeeping", length(Housekeeping_1[,4])),   rep("HeartSpecific", length(Specific_2[,7]))   ), 
                  sampleRank2=c(  "Housekeeping",   "HeartSpecific"  ), 
                  colours2=c(   "Housekeeping"="red2",    "HeartSpecific"="blue2"   ), 
                  path2=AllResults_g,   fileName2= paste("5-BoxViolin-all",   "H3",  sep = "_") ,  
                  title2="Genes",   xLab2="Type",    yLab2="H3",   
                  height2=4.0,   width2=2,   Ymin2=0, Ymax2=2.0)    ## width = 1 + 4*0.5, height=5cm 



## T test and Wilcoxon test  (unpaired).
MyHypothesisTest_1 <- function(vector1, vector2, file1) {
  sink(file=file1)
  print("######################## Apply continuity correction in the normal approximation for the p-value. ###############################################")
  
  print("##################################################################################################################################for comparing boxplot")
  wilcoxTest_2 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,  exact=FALSE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTest_2  )
  cat( "\n\nExact p-value:", wilcoxTest_2$p.value, "\n\n\n\n\n" ) 
  
  print("##################################################################################################################################")
  wilcoxTest_4 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,  exact=FALSE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTest_4  )
  cat( "\n\nExact p-value:", wilcoxTest_4$p.value, "\n\n\n\n\n" )
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  wilcoxTest_6 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,  exact=FALSE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTest_6  )
  cat( "\n\nExact p-value:", wilcoxTest_6$p.value, "\n\n\n\n\n\n" )
  
  
  print("######################## Don't apply continuity correction in the normal approximation for the p-value. ###############################################")
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  wilcoxTestB_2 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,  exact=FALSE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTestB_2  )
  cat( "\n\nExact p-value:", wilcoxTestB_2$p.value, "\n\n\n\n\n" ) 
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  wilcoxTestB_4 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,  exact=FALSE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTestB_4  )
  cat( "\n\nExact p-value:", wilcoxTestB_4$p.value, "\n\n\n\n\n" )
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  wilcoxTestB_6 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,  exact=FALSE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTestB_6  )
  cat( "\n\nExact p-value:", wilcoxTestB_6$p.value, "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" )
  
  
  
  
  
  print("######################## T-test, var.equal=FALSE. ###############################################")
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  tTest_2 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,   var.equal=FALSE,  conf.level=0.95)
  print( tTest_2  )
  cat( "\n\nExact p-value:", tTest_2$p.value, "\n\n\n\n\n" )
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  tTest_4 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,   var.equal=FALSE,  conf.level=0.95)
  print( tTest_4  )
  cat( "\n\nExact p-value:", tTest_4$p.value, "\n\n\n\n\n" )
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  tTest_6 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,   var.equal=FALSE,  conf.level=0.95)
  print( tTest_6  )
  cat( "\n\nExact p-value:", tTest_6$p.value, "\n\n\n\n\n\n" )
  
  
  print("######################## T-test, var.equal=TRUE. ##############################################################################################")
  
  print("##################################################################################################################################")
  tTestB_2 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=FALSE,   var.equal=TRUE,  conf.level=0.95)
  print( tTestB_2  )
  cat( "\n\nExact p-value:", tTestB_2$p.value, "\n\n\n\n\n" )
  print("##################################################################################################################################")
  
  print("##################################################################################################################################")
  tTestB_4 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=FALSE,   var.equal=TRUE,  conf.level=0.95)
  print( tTestB_4 )
  cat( "\n\nExact p-value:", tTestB_4$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  tTestB_6 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=FALSE,   var.equal=TRUE,  conf.level=0.95)
  print( tTestB_6  )
  cat( "\n\nExact p-value:", tTestB_6$p.value, "\n\n\n\n\n" )
  
  sink()
}

MyHypothesisTest_1(vector1=Housekeeping_1[,3],   vector2=Specific_2[,3],    file1=paste(AllResults_g,  "/1-BoxViolin-all",    ".TPM.txt",        sep = "")   )   
MyHypothesisTest_1(vector1=Housekeeping_1[,4],   vector2=Specific_2[,4],    file1=paste(AllResults_g,  "/2-BoxViolin-all",    ".NTR-1kb.txt",    sep = "")   )   
MyHypothesisTest_1(vector1=Housekeeping_1[,5],   vector2=Specific_2[,5],    file1=paste(AllResults_g,  "/3-BoxViolin-all",    ".NTR-500bp.txt",  sep = "")   )   
MyHypothesisTest_1(vector1=Housekeeping_1[,6],   vector2=Specific_2[,6],    file1=paste(AllResults_g,  "/4-BoxViolin-all",    ".NTR-max.txt",    sep = "")   )   
MyHypothesisTest_1(vector1=Housekeeping_1[,7],   vector2=Specific_2[,7],    file1=paste(AllResults_g,  "/5-BoxViolin-all",    ".H3.txt",         sep = "")   )   













Housekeeping_1[1:10,]
Specific_2[1:10,]

WT_oneClass <- ceiling( length(Housekeeping_1[,3])/4 )  ## have been sorted by TPM
WT_oneClass
WT_oneClass2 <- length(Housekeeping_1[,3]) - WT_oneClass * 3
WT_xAxis <-  c( rep("Lowest", WT_oneClass),   rep("Low", WT_oneClass),  rep("High", WT_oneClass) ,  rep("Highest", ceiling(WT_oneClass2 ) )  ) 
length(WT_xAxis)
length(Housekeeping_1[,3])

class1 <- c( (WT_oneClass*0+1): (WT_oneClass*1))
class2 <- c( (WT_oneClass*1+1): (WT_oneClass*2))
class3 <- c( (WT_oneClass*2+1): (WT_oneClass*3))
class4 <- c( (WT_oneClass*3+1): (WT_oneClass*4-1))


MyBoxViolinPlot_1(vector2= log2( c( Housekeeping_1[class1,3],   Housekeeping_1[class2,3],    Housekeeping_1[class3,3],   Housekeeping_1[class4,3]   ) + 1 ) ,   
                  sampleType2= WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("6-1-BoxViolin-4classes",   "TPM",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="TPM",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=12)    ## width = 1 + 4*0.5, height=5cm 


MyBoxViolinPlot_1(vector2= c( Housekeeping_1[class1,4],   Housekeeping_1[class2,4],    Housekeeping_1[class3,4],   Housekeeping_1[class4,4]   )  ,   
                  sampleType2= WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("6-2-BoxViolin-4classes",   "NTR-1kb",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="NTR",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=1)    ## width = 1 + 4*0.5, height=5cm 



MyBoxViolinPlot_1(vector2= c( Housekeeping_1[class1,5],   Housekeeping_1[class2,5],    Housekeeping_1[class3,5],   Housekeeping_1[class4,5]   )  ,   
                  sampleType2= WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("6-3-BoxViolin-4classes",   "NTR-500bp",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="NTR",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=1)    ## width = 1 + 4*0.5, height=5cm 




MyBoxViolinPlot_1(vector2= c( Housekeeping_1[class1,6],   Housekeeping_1[class2,6],    Housekeeping_1[class3,6],   Housekeeping_1[class4,6]   )  ,   
                  sampleType2= WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("6-4-BoxViolin-4classes",   "NTR-max",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="NTR",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=1)    ## width = 1 + 4*0.5, height=5cm 


MyBoxViolinPlot_1(vector2= c( Housekeeping_1[class1,7],   Housekeeping_1[class2,7],    Housekeeping_1[class3,7],   Housekeeping_1[class4,7]   )  ,   
                  sampleType2= WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("6-5-BoxViolin-4classes",   "H3",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="H3",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=2)    ## width = 1 + 4*0.5, height=5cm 























B_WT_oneClass <- ceiling( length(Specific_2[,3])/4 )
B_WT_oneClass
B_WT_oneClass2 <- length(Specific_2[,3]) - B_WT_oneClass * 3
B_WT_xAxis <-  c( rep("Lowest", B_WT_oneClass),   rep("Low", B_WT_oneClass),  rep("High", B_WT_oneClass) ,  rep("Highest", ceiling(B_WT_oneClass2 ) )  ) 
length(B_WT_xAxis)
length(Specific_2[,3])


B_class1 <- c( (B_WT_oneClass*0+1): (B_WT_oneClass*1))
B_class2 <- c( (B_WT_oneClass*1+1): (B_WT_oneClass*2))
B_class3 <- c( (B_WT_oneClass*2+1): (B_WT_oneClass*3))
B_class4 <- c( (B_WT_oneClass*3+1): (B_WT_oneClass*4-2))


MyBoxViolinPlot_1(vector2= log2( c( Specific_2[B_class1,3],   Specific_2[B_class2,3],    Specific_2[B_class3,3],   Specific_2[B_class4,3]   ) + 1 ) ,   
                  sampleType2= B_WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("7-1-BoxViolin-4classes",   "TPM",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="TPM",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=12)    ## width = 1 + 4*0.5, height=5cm 


MyBoxViolinPlot_1(vector2= c( Specific_2[B_class1,4],   Specific_2[B_class2,4],    Specific_2[B_class3,4],   Specific_2[B_class4,4]   )  ,   
                  sampleType2= B_WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("7-2-BoxViolin-4classes",   "NTR-1kb",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="NTR",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=1)    ## width = 1 + 4*0.5, height=5cm 



MyBoxViolinPlot_1(vector2= c( Specific_2[B_class1,5],   Specific_2[B_class2,5],    Specific_2[B_class3,5],   Specific_2[B_class4,5]   )  ,   
                  sampleType2= B_WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("7-3-BoxViolin-4classes",   "NTR-500bp",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="NTR",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=1)    ## width = 1 + 4*0.5, height=5cm 




MyBoxViolinPlot_1(vector2= c( Specific_2[B_class1,6],   Specific_2[B_class2,6],    Specific_2[B_class3,6],   Specific_2[B_class4,6]   )  ,   
                  sampleType2= B_WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("7-4-BoxViolin-4classes",   "NTR-max",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="NTR",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=1)    ## width = 1 + 4*0.5, height=5cm 



MyBoxViolinPlot_1(vector2= c( Specific_2[B_class1,7],   Specific_2[B_class2,7],    Specific_2[B_class3,7],   Specific_2[B_class4,7]   )  ,   
                  sampleType2= B_WT_xAxis, 
                  sampleRank2=c(  "Lowest",   "Low", "High",  "Highest"), 
                  colours2= "red2", 
                  path2=AllResults_g,   fileName2= paste("7-5-BoxViolin-4classes",   "H3",  sep = "_") ,  
                  title2="Genes",   xLab2="TPM",    yLab2="H3",   
                  height2=4.0,   width2=3,   Ymin2=0, Ymax2=2)    ## width = 1 + 4*0.5, height=5cm 










dim(Housekeeping_1)
dim(Specific_2)

Housekeeping_1[1:10,]
Specific_2[1:10,]

type2 <- rep("Housekeeping", nrow(Housekeeping_1) )
Housekeeping_1A <- cbind(Housekeeping_1, type2 )
type2 <- rep("HeartSpecific", nrow(Specific_2) )
Specific_2A     <- cbind(Specific_2, type2)
dim(Housekeeping_1A)
dim(Specific_2A)

Housekeeping_1A[1:10,]
Specific_2A[1:10,]

merged1 <- rbind(Housekeeping_1A, Specific_2A)
dim(merged1)
merged1[1:10,]
merged1[2300:2310,]




##################################################################################################################################################
index1 <-  order(merged1[,3])  ## sorted by TPM
length(index1)
index1[1:10]
Matrix_sort1  <- merged1[index1,]
dim(Matrix_sort1)
Matrix_sort1[1:20,]


numRows_oneClass <- floor( nrow(Matrix_sort1)/10 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 ) 
INDEX_6 <- seq(from = numRows_oneClass*5+1, to = numRows_oneClass*6, by =1 )
INDEX_7 <- seq(from = numRows_oneClass*6+1, to = numRows_oneClass*7, by =1 )
INDEX_8 <- seq(from = numRows_oneClass*7+1, to = numRows_oneClass*8, by =1 )
INDEX_9 <- seq(from = numRows_oneClass*8+1, to = numRows_oneClass*9, by =1 ) 
INDEX_10<- seq(from = numRows_oneClass*9+1, to = nrow(Matrix_sort1), by =1 )
nrow(Matrix_sort1) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4)  + length(INDEX_5)  + length(INDEX_6) 
                         + length(INDEX_7)  + length(INDEX_8)  + length(INDEX_9)  + length(INDEX_10) )



MyBoxViolinPlot_1(vector2= log2( c( Matrix_sort1[INDEX_1,3], Matrix_sort1[INDEX_2,3],  Matrix_sort1[INDEX_3,3],  Matrix_sort1[INDEX_4,3], Matrix_sort1[INDEX_5,3], 
                              Matrix_sort1[INDEX_6,3], Matrix_sort1[INDEX_7,3],  Matrix_sort1[INDEX_8,3],  Matrix_sort1[INDEX_9,3], Matrix_sort1[INDEX_10,3] ) + 1 ),   
                  sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("8-1A-sortByTPM",   "TPM",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="TPM",   
                  height2=4.0,   width2=6.0,   Ymin2=2.5, Ymax2=12)    ## width = 1 + 10*0.5, height=5cm 





MyGroupedBoxplot_1(vector2= log2( c( Matrix_sort1[INDEX_1,3], Matrix_sort1[INDEX_2,3],  Matrix_sort1[INDEX_3,3],  Matrix_sort1[INDEX_4,3], Matrix_sort1[INDEX_5,3], 
                                    Matrix_sort1[INDEX_6,3], Matrix_sort1[INDEX_7,3],  Matrix_sort1[INDEX_8,3],  Matrix_sort1[INDEX_9,3], Matrix_sort1[INDEX_10,3] ) + 1 ),   
                  sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  fill2=Matrix_sort1$type2,
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  path2=AllResults_g,   fileName2= paste("8-1B-sortByTPM",   "TPM",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="TPM",   
                  height2=4.0,   width2=10,   Ymin2=2.5, Ymax2=12)    ## width = 1 + 10*0.5, height=5cm 

 

summary( Matrix_sort1$TPM[Matrix_sort1$type2 == 'Housekeeping'] )
summary( Matrix_sort1$TPM[Matrix_sort1$type2 == 'HeartSpecific'] )

MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_1,3])[Matrix_sort1[INDEX_1,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_1,3])[Matrix_sort1[INDEX_1,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-1",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_2,3])[Matrix_sort1[INDEX_2,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_2,3])[Matrix_sort1[INDEX_2,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-2",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_3,3])[Matrix_sort1[INDEX_3,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_3,3])[Matrix_sort1[INDEX_3,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-3",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_4,3])[Matrix_sort1[INDEX_4,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_4,3])[Matrix_sort1[INDEX_4,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-4",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_5,3])[Matrix_sort1[INDEX_5,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_5,3])[Matrix_sort1[INDEX_5,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-5",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_6,3])[Matrix_sort1[INDEX_6,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_6,3])[Matrix_sort1[INDEX_6,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-6",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_7,3])[Matrix_sort1[INDEX_7,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_7,3])[Matrix_sort1[INDEX_7,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-7",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_8,3])[Matrix_sort1[INDEX_8,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_8,3])[Matrix_sort1[INDEX_8,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-8",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_9,3])[Matrix_sort1[INDEX_9,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_9,3])[Matrix_sort1[INDEX_9,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-9",    ".TPM.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_10,3])[Matrix_sort1[INDEX_10,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_10,3])[Matrix_sort1[INDEX_10,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-1B-sortByTPM-10",    ".TPM.txt",         sep = "")   )   










MyBoxViolinPlot_1(vector2= c( Matrix_sort1[INDEX_1,4], Matrix_sort1[INDEX_2,4],  Matrix_sort1[INDEX_3,4],  Matrix_sort1[INDEX_4,4], Matrix_sort1[INDEX_5,4], 
                              Matrix_sort1[INDEX_6,4], Matrix_sort1[INDEX_7,4],  Matrix_sort1[INDEX_8,4],  Matrix_sort1[INDEX_9,4], Matrix_sort1[INDEX_10,4] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("8-2A-sortByTPM",   "NTR-1kb",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="NTR-1kb",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort1[INDEX_1,4], Matrix_sort1[INDEX_2,4],  Matrix_sort1[INDEX_3,4],  Matrix_sort1[INDEX_4,4], Matrix_sort1[INDEX_5,4], 
                                Matrix_sort1[INDEX_6,4], Matrix_sort1[INDEX_7,4],  Matrix_sort1[INDEX_8,4],  Matrix_sort1[INDEX_9,4], Matrix_sort1[INDEX_10,4] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                  rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                  rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                   fill2=Matrix_sort1$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("8-2B-sortByTPM",   "NTR-1kb",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="NTR-1kb",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

Matrix_sort1[1:20,]


MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_1,4])[Matrix_sort1[INDEX_1,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_1,4])[Matrix_sort1[INDEX_1,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-1",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_2,4])[Matrix_sort1[INDEX_2,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_2,4])[Matrix_sort1[INDEX_2,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-2",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_3,4])[Matrix_sort1[INDEX_3,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_3,4])[Matrix_sort1[INDEX_3,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-3",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_4,4])[Matrix_sort1[INDEX_4,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_4,4])[Matrix_sort1[INDEX_4,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-4",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_5,4])[Matrix_sort1[INDEX_5,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_5,4])[Matrix_sort1[INDEX_5,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-5",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_6,4])[Matrix_sort1[INDEX_6,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_6,4])[Matrix_sort1[INDEX_6,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-6",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_7,4])[Matrix_sort1[INDEX_7,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_7,4])[Matrix_sort1[INDEX_7,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-7",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_8,4])[Matrix_sort1[INDEX_8,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_8,4])[Matrix_sort1[INDEX_8,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-8",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_9,4])[Matrix_sort1[INDEX_9,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_9,4])[Matrix_sort1[INDEX_9,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-9",    ".NTR-1kb.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_10,4])[Matrix_sort1[INDEX_10,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_10,4])[Matrix_sort1[INDEX_10,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-2B-sortByTPM-10",    ".NTR-1kb.txt",         sep = "")   )   











MyBoxViolinPlot_1(vector2= c( Matrix_sort1[INDEX_1,5], Matrix_sort1[INDEX_2,5],  Matrix_sort1[INDEX_3,5],  Matrix_sort1[INDEX_4,5], Matrix_sort1[INDEX_5,5], 
                              Matrix_sort1[INDEX_6,5], Matrix_sort1[INDEX_7,5],  Matrix_sort1[INDEX_8,5],  Matrix_sort1[INDEX_9,5], Matrix_sort1[INDEX_10,5] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("8-3A-sortByTPM",   "NTR-500bp",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="NTR-500bp",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort1[INDEX_1,5], Matrix_sort1[INDEX_2,5],  Matrix_sort1[INDEX_3,5],  Matrix_sort1[INDEX_4,5], Matrix_sort1[INDEX_5,5], 
                                Matrix_sort1[INDEX_6,5], Matrix_sort1[INDEX_7,5],  Matrix_sort1[INDEX_8,5],  Matrix_sort1[INDEX_9,5], Matrix_sort1[INDEX_10,5] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                  rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                  rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                   fill2=Matrix_sort1$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("8-3B-sortByTPM",   "NTR-500bp",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="NTR-500bp",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 




MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_1,5])[Matrix_sort1[INDEX_1,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_1,5])[Matrix_sort1[INDEX_1,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-1",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_2,5])[Matrix_sort1[INDEX_2,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_2,5])[Matrix_sort1[INDEX_2,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-2",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_3,5])[Matrix_sort1[INDEX_3,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_3,5])[Matrix_sort1[INDEX_3,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-3",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_4,5])[Matrix_sort1[INDEX_4,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_4,5])[Matrix_sort1[INDEX_4,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-4",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_5,5])[Matrix_sort1[INDEX_5,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_5,5])[Matrix_sort1[INDEX_5,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-5",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_6,5])[Matrix_sort1[INDEX_6,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_6,5])[Matrix_sort1[INDEX_6,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-6",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_7,5])[Matrix_sort1[INDEX_7,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_7,5])[Matrix_sort1[INDEX_7,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-7",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_8,5])[Matrix_sort1[INDEX_8,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_8,5])[Matrix_sort1[INDEX_8,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-8",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_9,5])[Matrix_sort1[INDEX_9,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_9,5])[Matrix_sort1[INDEX_9,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-9",    ".NTR-500bp.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_10,5])[Matrix_sort1[INDEX_10,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_10,5])[Matrix_sort1[INDEX_10,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-3B-sortByTPM-10",    ".NTR-500bp.txt",         sep = "")   )   












MyBoxViolinPlot_1(vector2= c( Matrix_sort1[INDEX_1,6], Matrix_sort1[INDEX_2,6],  Matrix_sort1[INDEX_3,6],  Matrix_sort1[INDEX_4,6], Matrix_sort1[INDEX_5,6], 
                              Matrix_sort1[INDEX_6,6], Matrix_sort1[INDEX_7,6],  Matrix_sort1[INDEX_8,6],  Matrix_sort1[INDEX_9,6], Matrix_sort1[INDEX_10,6] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("8-4A-sortByTPM",   "NTR-max",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="NTR-max",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort1[INDEX_1,6], Matrix_sort1[INDEX_2,6],  Matrix_sort1[INDEX_3,6],  Matrix_sort1[INDEX_4,6], Matrix_sort1[INDEX_5,6], 
                                Matrix_sort1[INDEX_6,6], Matrix_sort1[INDEX_7,6],  Matrix_sort1[INDEX_8,6],  Matrix_sort1[INDEX_9,6], Matrix_sort1[INDEX_10,6] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                  rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                  rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                   fill2=Matrix_sort1$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("8-4B-sortByTPM",   "NTR-max",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="NTR-max",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 





MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_1,6])[Matrix_sort1[INDEX_1,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_1,6])[Matrix_sort1[INDEX_1,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-1",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_2,6])[Matrix_sort1[INDEX_2,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_2,6])[Matrix_sort1[INDEX_2,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-2",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_3,6])[Matrix_sort1[INDEX_3,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_3,6])[Matrix_sort1[INDEX_3,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-3",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_4,6])[Matrix_sort1[INDEX_4,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_4,6])[Matrix_sort1[INDEX_4,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-4",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_5,6])[Matrix_sort1[INDEX_5,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_5,6])[Matrix_sort1[INDEX_5,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-5",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_6,6])[Matrix_sort1[INDEX_6,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_6,6])[Matrix_sort1[INDEX_6,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-6",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_7,6])[Matrix_sort1[INDEX_7,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_7,6])[Matrix_sort1[INDEX_7,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-7",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_8,6])[Matrix_sort1[INDEX_8,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_8,6])[Matrix_sort1[INDEX_8,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-8",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_9,6])[Matrix_sort1[INDEX_9,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_9,6])[Matrix_sort1[INDEX_9,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-9",    ".NTR-max.txt",         sep = "")   )   
MyHypothesisTest_1(vector1 = (Matrix_sort1[INDEX_10,6])[Matrix_sort1[INDEX_10,8] == 'Housekeeping'],   vector2 = (Matrix_sort1[INDEX_10,6])[Matrix_sort1[INDEX_10,8] == 'HeartSpecific'],    file1=paste(AllResults_g,  "/8-4B-sortByTPM-10",    ".NTR-max.txt",         sep = "")   )   






MyBoxViolinPlot_1(vector2= c( Matrix_sort1[INDEX_1,7], Matrix_sort1[INDEX_2,7],  Matrix_sort1[INDEX_3,7],  Matrix_sort1[INDEX_4,7], Matrix_sort1[INDEX_5,7], 
                              Matrix_sort1[INDEX_6,7], Matrix_sort1[INDEX_7,7],  Matrix_sort1[INDEX_8,7],  Matrix_sort1[INDEX_9,7], Matrix_sort1[INDEX_10,7] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("8-5A-sortByTPM",   "H3",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="H3",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=2)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort1[INDEX_1,7], Matrix_sort1[INDEX_2,7],  Matrix_sort1[INDEX_3,7],  Matrix_sort1[INDEX_4,7], Matrix_sort1[INDEX_5,7], 
                                Matrix_sort1[INDEX_6,7], Matrix_sort1[INDEX_7,7],  Matrix_sort1[INDEX_8,7],  Matrix_sort1[INDEX_9,7], Matrix_sort1[INDEX_10,7] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX_1)),     rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                  rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                  rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                   fill2=Matrix_sort1$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("8-5B-sortByTPM",   "H3",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="TPM (%)",    yLab2="H3",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=2)    ## width = 1 + 10*0.5, height=5cm 




















##################################################################################################################################################
index2 <-  order(merged1[,6])  ## sorted by NTR-max
length(index2)
index2[1:10]
Matrix_sort2  <- merged1[index2,]
dim(Matrix_sort2)
Matrix_sort2[1:20,]


numRows_oneClass <- floor( nrow(Matrix_sort2)/10 )
INDEX2_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX2_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX2_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX2_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX2_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 ) 
INDEX2_6 <- seq(from = numRows_oneClass*5+1, to = numRows_oneClass*6, by =1 )
INDEX2_7 <- seq(from = numRows_oneClass*6+1, to = numRows_oneClass*7, by =1 )
INDEX2_8 <- seq(from = numRows_oneClass*7+1, to = numRows_oneClass*8, by =1 )
INDEX2_9 <- seq(from = numRows_oneClass*8+1, to = numRows_oneClass*9, by =1 ) 
INDEX2_10<- seq(from = numRows_oneClass*9+1, to = nrow(Matrix_sort2), by =1 )
nrow(Matrix_sort2) - (length(INDEX2_1) + length(INDEX2_2) + length(INDEX2_3) + length(INDEX2_4)  + length(INDEX2_5)  + length(INDEX2_6) 
                      + length(INDEX2_7)  + length(INDEX2_8)  + length(INDEX2_9)  + length(INDEX2_10) )



MyBoxViolinPlot_1(vector2= log2( c( Matrix_sort2[INDEX2_1,3], Matrix_sort2[INDEX2_2,3],  Matrix_sort2[INDEX2_3,3],  Matrix_sort2[INDEX2_4,3], Matrix_sort2[INDEX2_5,3], 
                                    Matrix_sort2[INDEX2_6,3], Matrix_sort2[INDEX2_7,3],  Matrix_sort2[INDEX2_8,3],  Matrix_sort2[INDEX2_9,3], Matrix_sort2[INDEX2_10,3] ) + 1 ),   
                  sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                 rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                 rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("9-1A-sortByNTRmax",   "TPM",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="TPM",   
                  height2=4.0,   width2=6.0,   Ymin2=2.5, Ymax2=12)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2= log2( c( Matrix_sort2[INDEX2_1,3], Matrix_sort2[INDEX2_2,3],  Matrix_sort2[INDEX2_3,3],  Matrix_sort2[INDEX2_4,3], Matrix_sort2[INDEX2_5,3], 
                                     Matrix_sort2[INDEX2_6,3], Matrix_sort2[INDEX2_7,3],  Matrix_sort2[INDEX2_8,3],  Matrix_sort2[INDEX2_9,3], Matrix_sort2[INDEX2_10,3] ) + 1 ),   
                   sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                  rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                  rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                   fill2=Matrix_sort2$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("9-1B-sortByNTRmax",   "TPM",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="TPM",   
                   height2=4.0,   width2=10,   Ymin2=2.5, Ymax2=12)    ## width = 1 + 10*0.5, height=5cm 



summary( Matrix_sort2$TPM[Matrix_sort2$type2 == 'Housekeeping'] )
summary( Matrix_sort2$TPM[Matrix_sort2$type2 == 'HeartSpecific'] )





MyBoxViolinPlot_1(vector2= c( Matrix_sort2[INDEX2_1,4], Matrix_sort2[INDEX2_2,4],  Matrix_sort2[INDEX2_3,4],  Matrix_sort2[INDEX2_4,4], Matrix_sort2[INDEX2_5,4], 
                              Matrix_sort2[INDEX2_6,4], Matrix_sort2[INDEX2_7,4],  Matrix_sort2[INDEX2_8,4],  Matrix_sort2[INDEX2_9,4], Matrix_sort2[INDEX2_10,4] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                 rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                 rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("9-2A-sortByNTRmax",   "NTR-1kb",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="NTR-1kb",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort2[INDEX2_1,4], Matrix_sort2[INDEX2_2,4],  Matrix_sort2[INDEX2_3,4],  Matrix_sort2[INDEX2_4,4], Matrix_sort2[INDEX2_5,4], 
                                Matrix_sort2[INDEX2_6,4], Matrix_sort2[INDEX2_7,4],  Matrix_sort2[INDEX2_8,4],  Matrix_sort2[INDEX2_9,4], Matrix_sort2[INDEX2_10,4] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                  rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                  rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                   fill2=Matrix_sort2$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("9-2B-sortByNTRmax",   "NTR-1kb",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="NTR-1kb",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

Matrix_sort2[1:20,]






MyBoxViolinPlot_1(vector2= c( Matrix_sort2[INDEX2_1,5], Matrix_sort2[INDEX2_2,5],  Matrix_sort2[INDEX2_3,5],  Matrix_sort2[INDEX2_4,5], Matrix_sort2[INDEX2_5,5], 
                              Matrix_sort2[INDEX2_6,5], Matrix_sort2[INDEX2_7,5],  Matrix_sort2[INDEX2_8,5],  Matrix_sort2[INDEX2_9,5], Matrix_sort2[INDEX2_10,5] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                 rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                 rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("9-3A-sortByNTRmax",   "NTR-500bp",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="NTR-500bp",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort2[INDEX2_1,5], Matrix_sort2[INDEX2_2,5],  Matrix_sort2[INDEX2_3,5],  Matrix_sort2[INDEX2_4,5], Matrix_sort2[INDEX2_5,5], 
                                Matrix_sort2[INDEX2_6,5], Matrix_sort2[INDEX2_7,5],  Matrix_sort2[INDEX2_8,5],  Matrix_sort2[INDEX2_9,5], Matrix_sort2[INDEX2_10,5] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                  rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                  rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                   fill2=Matrix_sort2$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("9-3B-sortByNTRmax",   "NTR-500bp",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="NTR-1500bp",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 






MyBoxViolinPlot_1(vector2= c( Matrix_sort2[INDEX2_1,6], Matrix_sort2[INDEX2_2,6],  Matrix_sort2[INDEX2_3,6],  Matrix_sort2[INDEX2_4,6], Matrix_sort2[INDEX2_5,6], 
                              Matrix_sort2[INDEX2_6,6], Matrix_sort2[INDEX2_7,6],  Matrix_sort2[INDEX2_8,6],  Matrix_sort2[INDEX2_9,6], Matrix_sort2[INDEX2_10,6] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                 rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                 rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("9-4A-sortByNTRmax",   "NTR-max",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="NTR-max",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort2[INDEX2_1,6], Matrix_sort2[INDEX2_2,6],  Matrix_sort2[INDEX2_3,6],  Matrix_sort2[INDEX2_4,6], Matrix_sort2[INDEX2_5,6], 
                                Matrix_sort2[INDEX2_6,6], Matrix_sort2[INDEX2_7,6],  Matrix_sort2[INDEX2_8,6],  Matrix_sort2[INDEX2_9,6], Matrix_sort2[INDEX2_10,6] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                  rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                  rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                   fill2=Matrix_sort2$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("9-4B-sortByNTRmax",   "NTR-max",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="NTR-max",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=0.7)    ## width = 1 + 10*0.5, height=5cm 





MyBoxViolinPlot_1(vector2= c( Matrix_sort2[INDEX2_1,7], Matrix_sort2[INDEX2_2,7],  Matrix_sort2[INDEX2_3,7],  Matrix_sort2[INDEX2_4,7], Matrix_sort2[INDEX2_5,7], 
                              Matrix_sort2[INDEX2_6,7], Matrix_sort2[INDEX2_7,7],  Matrix_sort2[INDEX2_8,7],  Matrix_sort2[INDEX2_9,7], Matrix_sort2[INDEX2_10,7] ) ,   
                  sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                 rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                 rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red4",  "20-30"="orange",  "30-40"="orange4",  "40-50"="yellow", 
                               "50-60"="yellow4", "60-70"="green", "70-80"="green4", "80-90"="blue", "90-100"="blue4"), 
                  path2=AllResults_g,   fileName2= paste("9-5A-sortByNTRmax",   "H3",  sep = "_") , 
                  title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="H3",   
                  height2=4.0,   width2=6.0,   Ymin2=0, Ymax2=2)    ## width = 1 + 10*0.5, height=5cm 

MyGroupedBoxplot_1(vector2=  c( Matrix_sort2[INDEX2_1,7], Matrix_sort2[INDEX2_2,7],  Matrix_sort2[INDEX2_3,7],  Matrix_sort2[INDEX2_4,7], Matrix_sort2[INDEX2_5,7], 
                                Matrix_sort2[INDEX2_6,7], Matrix_sort2[INDEX2_7,7],  Matrix_sort2[INDEX2_8,7],  Matrix_sort2[INDEX2_9,7], Matrix_sort2[INDEX2_10,7] ) ,   
                   sampleType2=c( rep("0-10", length(INDEX2_1)),     rep("10-20", length(INDEX2_2)),  rep("20-30", length(INDEX2_3)),  rep("30-40", length(INDEX2_4)),
                                  rep("40-50", length(INDEX2_5)),    rep("50-60", length(INDEX2_6)),  rep("60-70", length(INDEX2_7)),  rep("70-80", length(INDEX2_8)),   
                                  rep("80-90", length(INDEX2_9)),    rep("90-100", length(INDEX2_10))   ), 
                   fill2=Matrix_sort2$type2,
                   sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", "50-60", "60-70", "70-80", "80-90", "90-100"),     
                   path2=AllResults_g,   fileName2= paste("9-5B-sortByNTRmax",   "H3",  sep = "_") , 
                   title2="housekeeping and heart-specific genes",   xLab2="NTR (%)",    yLab2="H3",   
                   height2=4.0,   width2=10,   Ymin2=0, Ymax2=2)    ## width = 1 + 10*0.5, height=5cm 






















