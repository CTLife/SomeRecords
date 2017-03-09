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
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 2) +   scale_x_discrete(limits=sampleRank2)  +
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "red", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-10adjust",          sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp3 <- ggplot(DataFrame_Local, aes(x=sampleType) ) + 
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray")   + scale_x_discrete(limits=sampleRank2)  +
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),    width = 0.1, size=0.3, colour = "black") +
    geom_boxplot( aes(y=yAxis),   width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = "red", notch=FALSE,  notchwidth = 0.15, alpha=1) + 
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
    geom_violin(aes(y=yAxis), fill = "gray", colour = "gray", adjust = 10) + scale_x_discrete(limits=sampleRank2)  +  
    #geom_errorbar(aes(ymin=min,ymax=max),  data=whisk_1(DataFrame_Local),   width = 0.1, size=0.1, colour = "black") +
    geom_boxplot( aes(y=yAxis),  width=0.3, size=0.5, fill=NA, outlier.size=0,  colour = colours2, notch=FALSE,  notchwidth = 0.15, alpha=1) + 
    stat_summary( aes(y=yAxis),   fun.y=mean, colour="white", geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1(textSize1=14, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1(ggplot2Figure1=FigureTemp5,  path1=path2, fileName1=paste(fileName2, "_ViolinPlot-10Adjust-colour",   sep="",  collapse=NULL),  height1=height2, width1=width2)
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



rawMatrix1  <- read.table("5-NTR-TPM/Genes_NTR1kb_NTR500bp_TPM",  header=TRUE,   sep="",   quote = "",   comment.char = "")  
dim(rawMatrix1)

TPM1 <- rawMatrix1[, 5]
NTR1 <- rawMatrix1[, 13]
H3_1 <- rawMatrix1[, 15]
length(TPM1)
length(NTR1) 
length(H3_1) 







index1 <-  order(H3_1)  ## sorted by H3
length(index1)
index1[1:10]
Matrix_sort_H3  <- rawMatrix1[index1,]
dim(Matrix_sort_H3)
H3_1[index1[1:10]]

numRows_oneClass <- floor( nrow(Matrix_sort_H3)/10 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 ) 
INDEX_6 <- seq(from = numRows_oneClass*5+1, to = numRows_oneClass*6, by =1 )
INDEX_7 <- seq(from = numRows_oneClass*6+1, to = numRows_oneClass*7, by =1 )
INDEX_8 <- seq(from = numRows_oneClass*7+1, to = numRows_oneClass*8, by =1 )
INDEX_9 <- seq(from = numRows_oneClass*8+1, to = numRows_oneClass*9, by =1 ) 
INDEX_10<- seq(from = numRows_oneClass*9+1, to = numRows_oneClass*10,by =1 )
nrow(Matrix_sort_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4)  + length(INDEX_5)  + length(INDEX_6) 
                         + length(INDEX_7)  + length(INDEX_8)  + length(INDEX_9)  + length(INDEX_10) )

TPM2 <- Matrix_sort_H3[, 5]

MyBoxViolinPlot_1(vector2= log2( c( TPM2[INDEX_1], TPM2[INDEX_2],  TPM2[INDEX_3],  TPM2[INDEX_4], TPM2[INDEX_5], 
                              TPM2[INDEX_6], TPM2[INDEX_7],  TPM2[INDEX_8],  TPM2[INDEX_9], TPM2[INDEX_10] ) + 1 ),   
                  sampleType2=c( rep("0-10", length(INDEX_1)),    rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", 
                                 "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red",  "20-30"="red",  "30-40"="red",  "40-50"="red", 
                               "50-60"="red", "60-70"="red", "70-80"="red", "80-90"="red", "90-100"="red"), 
                  path2="z-figures",   fileName2= paste("1-TPM-BoxViolin-byH3",   "500bp",  sep = "_") ,  
                  title2="H3-TPM",   xLab2="H3 (%)",    yLab2="log2(TPM+1)",   
                  height2=4.0,   width2=5.0,   Ymin2=0, Ymax2=12.3)    ## width = 1 + 10*0.5, height=5cm 





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

## T test and Wilcoxon test  (paired)
MyHypothesisTest_2 <- function(vector1, vector2, file1) {
  sink(file=file1)
  
  print("######################## Apply continuity correction in the normal approximation for the p-value. ###############################################")
  print("##################################################################################################################################")
  wilcoxTest_1 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,   exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTest_1  )
  cat( "\n\nExact p-value:", wilcoxTest_1$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  wilcoxTest_3 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,   exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTest_3  )
  cat( "\n\nExact p-value:", wilcoxTest_3$p.value, "\n\n\n\n\n" ) 
  
  print("##################################################################################################################################")
  wilcoxTest_5 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,   exact=TRUE, correct=TRUE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTest_5  )
  cat( "\n\nExact p-value:", wilcoxTest_5$p.value, "\n\n\n\n\n" )
  
  
  print("######################## Don't apply continuity correction in the normal approximation for the p-value. ###############################################")
  print("##################################################################################################################################")
  wilcoxTestB_1 <- wilcox.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,   exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTestB_1  )
  cat( "\n\nExact p-value:", wilcoxTestB_1$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  wilcoxTestB_3 <- wilcox.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,   exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTestB_3 )
  cat( "\n\nExact p-value:", wilcoxTestB_3$p.value, "\n\n\n\n\n" ) 
  
  print("##################################################################################################################################")
  wilcoxTestB_5 <- wilcox.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,   exact=TRUE, correct=FALSE,  conf.int=TRUE,  conf.level=0.95)
  print( wilcoxTestB_5  )
  cat( "\n\nExact p-value:", wilcoxTestB_5$p.value, "\n\n\n\n\n" )
  
  
  
  
  print("######################## T-test, var.equal=FALSE. ###############################################")
  print("##################################################################################################################################")
  tTest_1 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,    var.equal=FALSE,  conf.level=0.95)
  print( tTest_1  )
  cat( "\n\nExact p-value:", tTest_1$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  tTest_3 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,    var.equal=FALSE,  conf.level=0.95)
  print( tTest_3  )
  cat( "\n\nExact p-value:", tTest_3$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  tTest_5 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,    var.equal=FALSE,  conf.level=0.95)
  print( tTest_5  )
  cat( "\n\nExact p-value:", tTest_5$p.value, "\n\n\n\n\n" )
  
  
  print("######################## T-test, var.equal=TRUE. ##############################################################################################")
  print("##################################################################################################################################")
  tTestB_1 <- t.test(x=vector1, y=vector2, alternative="two.sided",  mu=0,   paired=TRUE,    var.equal=TRUE,  conf.level=0.95)
  print( tTestB_1  )
  cat( "\n\nExact p-value:", tTestB_1$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  tTestB_3 <- t.test(x=vector1, y=vector2, alternative="less",       mu=0,   paired=TRUE,    var.equal=TRUE,  conf.level=0.95)
  print( tTestB_3  )
  cat( "\n\nExact p-value:", tTestB_3$p.value, "\n\n\n\n\n" )
  
  print("##################################################################################################################################")
  tTestB_5 <- t.test(x=vector1, y=vector2, alternative="greater",    mu=0,   paired=TRUE,    var.equal=TRUE,  conf.level=0.95)
  print( tTestB_5  )
  cat( "\n\nExact p-value:", tTestB_5$p.value, "\n\n\n\n\n" )
  
  sink() 
}




MyHypothesisTest_1(vector1=TPM2[INDEX_1],  vector2=TPM2[INDEX_2],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-1-vs-2.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_2],  vector2=TPM2[INDEX_3],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-2-vs-3.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_3],  vector2=TPM2[INDEX_4],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-3-vs-4.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_4],  vector2=TPM2[INDEX_5],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-4-vs-5.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_5],  vector2=TPM2[INDEX_6],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-5-vs-6.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_6],  vector2=TPM2[INDEX_7],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-6-vs-7.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_7],  vector2=TPM2[INDEX_8],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-7-vs-8.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_8],  vector2=TPM2[INDEX_9],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-8-vs-9.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_9],  vector2=TPM2[INDEX_10],  file1=paste("z-figures",  "/1-TPM-BoxViolin-byH3-9-vs-10.unpaired.txt",     sep = "")   )   












index1 <-  runif(24412, 1, 24412 )  ## random order
length(index1)
index1[1:10]
Matrix_sort_H3  <- rawMatrix1[index1,]
dim(Matrix_sort_H3)
H3_1[index1[1:10]]

numRows_oneClass <- floor( nrow(Matrix_sort_H3)/10 )
INDEX_1 <- seq(from = numRows_oneClass*0+1, to = numRows_oneClass*1, by =1 ) 
INDEX_2 <- seq(from = numRows_oneClass*1+1, to = numRows_oneClass*2, by =1 )
INDEX_3 <- seq(from = numRows_oneClass*2+1, to = numRows_oneClass*3, by =1 )
INDEX_4 <- seq(from = numRows_oneClass*3+1, to = numRows_oneClass*4, by =1 )
INDEX_5 <- seq(from = numRows_oneClass*4+1, to = numRows_oneClass*5, by =1 ) 
INDEX_6 <- seq(from = numRows_oneClass*5+1, to = numRows_oneClass*6, by =1 )
INDEX_7 <- seq(from = numRows_oneClass*6+1, to = numRows_oneClass*7, by =1 )
INDEX_8 <- seq(from = numRows_oneClass*7+1, to = numRows_oneClass*8, by =1 )
INDEX_9 <- seq(from = numRows_oneClass*8+1, to = numRows_oneClass*9, by =1 ) 
INDEX_10<- seq(from = numRows_oneClass*9+1, to = numRows_oneClass*10,by =1 )
nrow(Matrix_sort_H3) - (length(INDEX_1) + length(INDEX_2) + length(INDEX_3) + length(INDEX_4)  + length(INDEX_5)  + length(INDEX_6) 
                        + length(INDEX_7)  + length(INDEX_8)  + length(INDEX_9)  + length(INDEX_10) )

TPM2 <- Matrix_sort_H3[, 5]

MyBoxViolinPlot_1(vector2= log2( c( TPM2[INDEX_1], TPM2[INDEX_2],  TPM2[INDEX_3],  TPM2[INDEX_4], TPM2[INDEX_5], 
                                    TPM2[INDEX_6], TPM2[INDEX_7],  TPM2[INDEX_8],  TPM2[INDEX_9], TPM2[INDEX_10] ) + 1 ),   
                  sampleType2=c( rep("0-10", length(INDEX_1)),    rep("10-20", length(INDEX_2)),  rep("20-30", length(INDEX_3)),  rep("30-40", length(INDEX_4)),
                                 rep("40-50", length(INDEX_5)),    rep("50-60", length(INDEX_6)),  rep("60-70", length(INDEX_7)),  rep("70-80", length(INDEX_8)),   
                                 rep("80-90", length(INDEX_9)),    rep("90-100", length(INDEX_10))   ), 
                  sampleRank2=c( "0-10",  "10-20",  "20-30",  "30-40",  "40-50", 
                                 "50-60", "60-70", "70-80", "80-90", "90-100"),     
                  colours2=c(  "0-10"="red",  "10-20"="red",  "20-30"="red",  "30-40"="red",  "40-50"="red", 
                               "50-60"="red", "60-70"="red", "70-80"="red", "80-90"="red", "90-100"="red"), 
                  path2="z-figures",   fileName2= paste("2-TPM-BoxViolin-byRandom",   "500bp",  sep = "_") ,  
                  title2="H3-TPM",   xLab2="H3 (%)",    yLab2="log2(TPM+1)",   
                  height2=4.0,   width2=5.0,   Ymin2=0, Ymax2=12.3)    ## width = 1 + 10*0.5, height=5cm 






MyHypothesisTest_1(vector1=TPM2[INDEX_1],  vector2=TPM2[INDEX_2],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-1-vs-2.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_2],  vector2=TPM2[INDEX_3],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-2-vs-3.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_3],  vector2=TPM2[INDEX_4],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-3-vs-4.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_4],  vector2=TPM2[INDEX_5],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-4-vs-5.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_5],  vector2=TPM2[INDEX_6],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-5-vs-6.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_6],  vector2=TPM2[INDEX_7],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-6-vs-7.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_7],  vector2=TPM2[INDEX_8],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-7-vs-8.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_8],  vector2=TPM2[INDEX_9],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-8-vs-9.unpaired.txt",     sep = "")   )   
MyHypothesisTest_1(vector1=TPM2[INDEX_9],  vector2=TPM2[INDEX_10],  file1=paste("z-figures",  "/2-TPM-BoxViolin-byRandom-9-vs-10.unpaired.txt",     sep = "")   )   



