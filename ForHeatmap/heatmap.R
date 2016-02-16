library("reshape2")
library("ggplot2") 
library("grid")
library("Cairo")
library("RColorBrewer")
library("gplots")  
library("stats")
library("KernSmooth")
library("matrixStats")
library("extrafont")
font_import()
fonttable()
loadfonts()
loadfonts(device="postscript")
names(postscriptFonts())







pTheme_1 <- function( textSize=12 ) {  
  theme(  
    line  = element_line(colour="black", size=1.0,   linetype=1,      lineend=NULL),                                                          ## all line elements.  局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black", size=1.0,   linetype=1,      fill="transparent" ),                                                   ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face=NULL,  colour="black",  size=textSize, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),   ## all text elements.  "serif" for a serif-serif font. 所有文本相关属性.
    title = element_text(family="serif",  face=NULL,  colour="black",  size=textSize, hjust=NULL, vjust=2,   angle=NULL, lineheight=NULL),    ## all title elements: plot, axes, legends. hjust:水平对齐的方向.  所有标题属性.
    
    axis.title        = element_text(family="serif", face=NULL, colour="black", size=textSize,    hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x      = element_text(family="serif", face=NULL, colour="black", size=textSize,    hjust=NULL, vjust=-0.5,  angle=NULL, lineheight=NULL),      ## x axis label (element_text; inherits from axis.title)
    axis.title.y      = element_text(family="serif", face=NULL, colour="black", size=textSize,    hjust=NULL, vjust=1.5,   angle=NULL, lineheight=NULL),      ## y axis label (element_text; inherits from axis.title)
    axis.text         = element_text(family="serif", face=NULL, colour="black", size=textSize-2,  hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x       = element_text(family="serif", face=NULL, colour="black", size=textSize-2,  hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y       = element_text(family="serif", face=NULL, colour="black", size=textSize-2,  hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    axis.ticks        = element_line(colour="black", size=0.3, linetype=1, lineend=NULL),                                                                     ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.3, linetype=1, lineend=NULL),                                                                     ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.3, linetype=1, lineend=NULL),                                                                     ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(1.5, "mm", data=NULL),                                                                                                     ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.ticks.margin = grid::unit(1.0, "mm", data=NULL),  	                                                                                                  ## space between tick mark and tick label (unit),  ‘"mm"’ Millimetres.  10 mm = 1 cm. 刻度线和刻度标签之间的间距.                                                                           
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                            ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                            ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	                                                              ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.margin        = grid::unit(0.3, "mm", data=NULL), 	                                                  ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(5, "mm", data=NULL) , 	                                                  ## size of legend keys (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(5.7, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(5, "mm", data=NULL) ,                                                     ## key background width (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize-2, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.3, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.margin       = grid::unit(2, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                              ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	  ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(8, 8, 8, 8), "mm", data=NULL), 	                                                                                  ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL, size=NULL, linetype=NULL, fill=NULL ), 	                                                        ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 










pTheme_2 <- function( textSize=12 ) {  
  theme(  
    line  = element_line(colour="black", size=1.0,   linetype=1,      lineend=NULL),                                                           ## all line elements.  局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black", size=1.0,   linetype=1,      fill="transparent" ),                                                    ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face=NULL,  colour="black",  size=textSize, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),    ## all text elements.  "serif" for a serif-serif font. 所有文本相关属性.
    title = element_text(family="serif",  face=NULL,  colour="black",  size=textSize, hjust=NULL, vjust=2,   angle=NULL, lineheight=NULL),     ## all title elements: plot, axes, legends. hjust:水平对齐的方向.  所有标题属性.
    
    axis.title        = element_text(family="serif", face=NULL, colour="black", size=textSize,    hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x      = element_text(family="serif", face=NULL, colour="black", size=textSize,    hjust=NULL, vjust=-0.5,  angle=NULL, lineheight=NULL),      ## x axis label (element_text; inherits from axis.title)
    axis.title.y      = element_text(family="serif", face=NULL, colour="black", size=textSize,    hjust=NULL, vjust=1.5,   angle=NULL, lineheight=NULL),      ## y axis label (element_text; inherits from axis.title)
    axis.text         = element_text(family="serif", face=NULL, colour="black", size=textSize-2,  hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x       = element_text(family="serif", face=NULL, colour="black", size=textSize-2,  hjust=1, vjust=1, angle=45, lineheight=NULL),               ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y       = element_text(family="serif", face=NULL, colour="black", size=textSize-2,  hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    axis.ticks        = element_line(colour="black", size=0.3, linetype=1, lineend=NULL),                                                                     ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.3, linetype=1, lineend=NULL),                                                                     ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.3, linetype=1, lineend=NULL),                                                                     ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(1.5, "mm", data=NULL),                                                                                                     ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.ticks.margin = grid::unit(1.0, "mm", data=NULL),                                                                                                     ## space between tick mark and tick label (unit),  ‘"mm"’ Millimetres.  10 mm = 1 cm. 刻度线和刻度标签之间的间距.                                                                           
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                            ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	                                                            ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.margin        = grid::unit(0.5, "mm", data=NULL), 	                                                  ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(5, "mm", data=NULL) , 	                                                  ## size of legend keys (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(5.7, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(5, "mm", data=NULL) ,                                                     ## key background width (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize-2, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	          ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	    ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	              ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	              ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.3, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.margin       = grid::unit(2, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                              ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	  ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(8, 8, 8, 8), "mm", data=NULL), 	                                                                                  ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL, size=NULL, linetype=NULL, fill=NULL ), 	                                                        ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 








normalize_1 <- function(vector1) {
  max1 = max(vector1)  
  min1 = min(vector1)  
  lower1 = 0  
  upper1 = 1
  vector2 = lower1 + (upper1-lower1)*(vector1-min1)/(max1-min1)
  return(vector2)
}









####################################################################  Start  ##########################################################################################################################################
matrix1 <- read.table("sorted_1097TFs",         header=TRUE,   sep="\t",    quote = "",   comment.char = "")
dim(matrix1)


matrix1A <- matrix1[, c(2:6)]
dim(matrix1A)
colnames(matrix1A) 
rownames(matrix1A) 

colnames(matrix1A) <- c("E9.5.Atria" ,   "E9.5.Outflowtrack",  "E9.5.Ventricle" ,   "E12.5" ,   "Adult"  ) 
rownames(matrix1A) <- paste( matrix1[,1] )
matrix1A
colnames(matrix1A) 
rownames(matrix1A) 


###################################
heatmap_1 <- as.matrix( matrix1A )   
rownames(heatmap_1)
colnames(heatmap_1)
dim(heatmap_1)
min(heatmap_1)
max(heatmap_1)

heatmap_1[1,]
normalize_1(heatmap_1[1,] )

heatmap_2  <- t( apply(heatmap_1, 1, normalize_1) )
dim(heatmap_2)
heatmap_2[1:5, ]
min(heatmap_2)
max(heatmap_2)
rownames(heatmap_2)
colnames(heatmap_2)
colnames(heatmap_2) <- c(  "E9.5_Ventricle", "E9.5_Atria", "E9.5_Outflowtract", "E12.5",  "Adult"  ) 

heatmap_3 <- melt( heatmap_2) 
rownames(heatmap_3)
colnames(heatmap_3)



## heatmap with lables
zp1 <- ggplot( heatmap_3, aes(x = Var2, y = Var1, fill = value) )
zp1 <- zp1 + geom_tile()  + xlab("Samples") + ylab("Genes") +  ggtitle("Expression Level of 1097 TFs") 
zp1 <- zp1 + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green4", midpoint = 0.4,  limits=c(0, 1) ) 
zp1 <- zp1 + pTheme_2(textSize=11) 

file1= "1-heatmap-Labels-yellowGreen"
postscript(file=paste(file1, ".eps", sep="", collapse=NULL),   height = 4.5,  width = 3.5,   family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)
print(zp1)
dev.off() 
ggsave(filename = paste(file1, ".svg", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".pdf", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".png", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )


## heatmap with lables
zp2 <- ggplot( heatmap_3, aes(x = Var2, y = Var1, fill = value) )
zp2 <- zp2 + geom_tile()  + xlab("Samples") + ylab("Genes") +  ggtitle("Expression Level of 1097 TFs") 
zp2 <- zp2 + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0.4,  limits=c(0, 1) ) 
zp2 <- zp2 + pTheme_2(textSize=11) 

file1= "2-heatmap-Labels-blueRed"
postscript(file=paste(file1, ".eps", sep="", collapse=NULL),   height = 4.5,  width = 3.5,   family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)
print(zp2)
dev.off() 
ggsave(filename = paste(file1, ".svg", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".pdf", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".png", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )


## heatmap with lables
zp3 <- ggplot( heatmap_3, aes(x = Var2, y = Var1, fill = value) )
zp3 <- zp3 + geom_tile()  + xlab("Samples") + ylab("Genes") +  ggtitle("Expression Level of 1097 TFs") 
zp3 <- zp3 + scale_fill_gradient2( low = "white",  high = "red", midpoint = 0.4,  limits=c(0, 1) ) 
zp3 <- zp3 + pTheme_2(textSize=11) 

file1= "3-heatmap-Labels-whiteRed"
postscript(file=paste(file1, ".eps", sep="", collapse=NULL),   height = 4.5,  width = 3.5,   family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)
print(zp3)
dev.off() 
ggsave(filename = paste(file1, ".svg", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".pdf", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".png", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )









## heatmap without lables
zp1 <- ggplot( heatmap_3, aes(x = Var2, y = as.numeric(Var1), fill = value) )    ## as.numeric() is nessesary.
zp1 <- zp1 + geom_tile()  + xlab("Samples") + ylab("Genes") +  ggtitle("Expression Level of 1097 TFs") 
zp1 <- zp1 + scale_fill_gradient2( low = "yellow", mid = "yellowgreen", high = "green4", midpoint = 0.4,  limits=c(0, 1) ) 
zp1 <- zp1 +  scale_x_discrete(expand = c(0, 0))  +  scale_y_continuous(expand = c(0, 0))  +  pTheme_2(textSize=11) 

file1= "1-heatmap-NoLabels-yellowGreen"
postscript(file=paste(file1, ".eps", sep="", collapse=NULL),   height = 4.5,  width = 3.5,   family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)
print(zp1)
dev.off() 
ggsave(filename = paste(file1, ".svg", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".pdf", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".png", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )


## heatmap without lables
zp2 <- ggplot( heatmap_3, aes(x = Var2, y = as.numeric(Var1), fill = value) )    ## as.numeric() is nessesary.
zp2 <- zp2 + geom_tile()  + xlab("Samples") + ylab("Genes") +  ggtitle("Expression Level of 1097 TFs") 
zp2 <- zp2 + scale_fill_gradient2( low = "blue", mid = "white", high = "red", midpoint = 0.4,  limits=c(0, 1) ) 
zp2 <- zp2 + scale_x_discrete(expand = c(0, 0))  +  scale_y_continuous(expand = c(0, 0))  + pTheme_2(textSize=11) 

file1= "2-heatmap-NoLabels-blueRed"
postscript(file=paste(file1, ".eps", sep="", collapse=NULL),   height = 4.5,  width = 3.5,   family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)
print(zp2)
dev.off() 
ggsave(filename = paste(file1, ".svg", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".pdf", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".png", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )


## heatmap without lables
zp3 <- ggplot( heatmap_3, aes(x = Var2, y = as.numeric(Var1), fill = value) )    ## as.numeric() is nessesary.
zp3 <- zp3 + geom_tile()  + xlab("Samples") + ylab("Genes") +  ggtitle("Expression Level of 1097 TFs") 
zp3 <- zp3 + scale_fill_gradient2( low = "white",  high = "red", midpoint = 0.4,  limits=c(0, 1) ) 
zp3 <- zp3 + scale_x_discrete(expand = c(0, 0))  +  scale_y_continuous(expand = c(0, 0))  + pTheme_2(textSize=11) 

file1= "3-heatmap-NoLabels-whiteRed"
postscript(file=paste(file1, ".eps", sep="", collapse=NULL),   height = 4.5,  width = 3.5,   family = "serif",  paper = "special",  onefile = FALSE,  horizontal = FALSE)
print(zp3)
dev.off() 
ggsave(filename = paste(file1, ".svg", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".pdf", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )
ggsave(filename = paste(file1, ".png", sep="", collapse=NULL),   height=4.5,  width=3.5,   dpi = 300 )





