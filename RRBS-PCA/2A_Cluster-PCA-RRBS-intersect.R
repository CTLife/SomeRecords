###########################################################################
FileLists <- list(
"1-Coverage-CpG/9_W1365C-boy-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/9_Q5-W1365D-boy-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/9_W1365F-father-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/9_W1365F-father-IVF-frozen_Rep2.bismark.cov",
"1-Coverage-CpG/9_Q1-W1365M-mother-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/10_W1733C-boy-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/10_W1733D-boy-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/10_W1733F-father-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/10_W1733M-Mother-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/11_W1398C-boy-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/11_W1398D-boy-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/11_W1398F-father-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/11_W1398M-mother-IVF-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/12_W1579C-boy-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/12_W1579C-boy-ICSI-fresh_Rep2.bismark.cov",
"1-Coverage-CpG/12_Q17-W1579D-boy-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/12_Q18-W1579F-father-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/12_Q19-W1579M-mother-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/13_W1647C-boy-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/13_W1647D-boy-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/13_W1647F-father-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/13_Q23-W1647M-mother-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/14_W1719C-boy-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/14_W1719D-boy-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/14_W1719F-father-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/14_Q15-W1719M-mother-ICSI-fresh_Rep1.bismark.cov",
"1-Coverage-CpG/15_Q4-W871D-boy-ICSI-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/15_Q6-W871C-boy-ICSI-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/15_Q7-W871F-father-ICSI-frozen_Rep1.bismark.cov",
"1-Coverage-CpG/15_Q10-W871M-mother-ICSI-frozen_Rep1.bismark.cov" 
)

mySampleID <- list("W1365-C-boy-IVF-frozen_1",         "W1365-D-boy-IVF-frozen_2", 
                   "W1365-F-father-IVF-frozen_3",      "W1365-F-father-IVF-frozen_4",   "W1365-M-mother-IVF-frozen_5",
                   "W1733-C-boy-IVF-frozen_6",         "W1733-D-boy-IVF-frozen_7", 
                   "W1733-F-father-IVF-frozen_8",      "W1733-M-mother-IVF-frozen_9",
                   "W1398-C-boy-IVF-frozen_10",         "W1398-D-boy-IVF-frozen_11", 
                   "W1398-F-father-IVF-frozen_12",      "W1398-M-mother-IVF-frozen_13",
                   "W1579-C-boy-ICSI-fresh_14",         "W1579-C-boy-ICSI-fresh_15",      "W1579-D-boy-ICSI-fresh_16", 
                   "W1579-F-father-ICSI-fresh_17",      "W1579-M-mother-ICSI-fresh_18",
                   "W1647-C-boy-ICSI-fresh_19",         "W1647-D-boy-ICSI-fresh_20", 
                   "W1647-F-father-ICSI-fresh_21",      "W1647-M-mother-ICSI-fresh_22",
                   "W1719-C-boy-ICSI-fresh_23",         "W1719-D-boy-ICSI-fresh_24", 
                   "W1719-F-father-ICSI-fresh_25",      "W1719-M-mother-ICSI-fresh_26",
                   "W871-C-boy-ICSI-frozen_27",         "W871-D-boy-ICSI-frozen_28", 
                   "W871-F-father-ICSI-frozen_29",      "W871-M-mother-ICSI-frozen_30" 
)

myOutDir <- "2A-all-CDFM-intersection"
myMinReads <- 5




myType1 <- c("W1365-C-boy-IVF-frozen_1"="red",    "W1365-D-boy-IVF-frozen_2"="red", 
             "W1365-F-father-IVF-frozen_3"="blue",       "W1365-F-father-IVF-frozen_4"="blue",   "W1365-M-mother-IVF-frozen_5"="cyan",
             "W1733-C-boy-IVF-frozen_6"="red",          "W1733-D-boy-IVF-frozen_7"="red", 
             "W1733-F-father-IVF-frozen_8"="blue",       "W1733-M-mother-IVF-frozen_9"="cyan",
             "W1398-C-boy-IVF-frozen_10"="red",         "W1398-D-boy-IVF-frozen_11"="red", 
             "W1398-F-father-IVF-frozen_12"="blue",      "W1398-M-mother-IVF-frozen_13"="cyan",
             "W1579-C-boy-ICSI-fresh_14"="red",         "W1579-C-boy-ICSI-fresh_15"="red",      "W1579-D-boy-ICSI-fresh_16"="red", 
             "W1579-F-father-ICSI-fresh_17"="blue",      "W1579-M-mother-ICSI-fresh_18"="cyan",
             "W1647-C-boy-ICSI-fresh_19"="red",         "W1647-D-boy-ICSI-fresh_20"="red", 
             "W1647-F-father-ICSI-fresh_21"="blue",      "W1647-M-mother-ICSI-fresh_22"="cyan",
             "W1719-C-boy-ICSI-fresh_23"="red",         "W1719-D-boy-ICSI-fresh_24"="red", 
             "W1719-F-father-ICSI-fresh_25"="blue",      "W1719-M-mother-ICSI-fresh_26"="cyan",
             "W871-C-boy-ICSI-frozen_27"="red",         "W871-D-boy-ICSI-frozen_28"="red", 
             "W871-F-father-ICSI-frozen_29"="blue",      "W871-M-mother-ICSI-frozen_30"="cyan" 
)



myType2 <- c("boy" ,        "boy" , 
             "father" ,     "father" ,      "mother" , 
             "boy" ,        "boy" , 
             "father" ,     "mother" ,  
             "boy" ,        "boy" , 
             "father" ,     "mother" , 
             "boy" ,        "boy" ,       "boy" , 
             "father" ,     "mother" , 
             "boy" ,        "boy" , 
             "father" ,     "mother" , 
             "boy" ,        "boy" , 
             "father" ,     "mother" , 
             "boy" ,        "boy" , 
             "father" ,     "mother" 
 
  )


myType3 <- c("IVF-frozen" ,     "IVF-frozen" , 
             "IVF-frozen" ,     "IVF-frozen" ,   "IVF-frozen" , 
             "IVF-frozen" ,     "IVF-frozen" , 
             "IVF-frozen" ,     "IVF-frozen" , 
             "IVF-frozen" ,     "IVF-frozen" , 
             "IVF-frozen" ,     "IVF-frozen" , 
             "ICSI-fresh" ,     "ICSI-fresh",  "ICSI-fresh", 
             "ICSI-fresh" ,     "ICSI-fresh", 
             "ICSI-fresh" ,     "ICSI-fresh", 
             "ICSI-fresh" ,     "ICSI-fresh", 
             "ICSI-fresh" ,     "ICSI-fresh",
             "ICSI-fresh" ,     "ICSI-fresh", 
             "ICSI-frozen" ,     "ICSI-frozen",
             "ICSI-frozen" ,     "ICSI-frozen"  
)




myType2_shape <- c("boy"=19 ,        "boy"=19 , 
             "father" =4,     "father" =4,      "mother"=0 , 
             "boy" =19,        "boy"=19 , 
             "father"=4 ,     "mother"=0 ,  
             "boy"=19 ,        "boy" =19, 
             "father"=4 ,     "mother"=0 , 
             "boy" =19,        "boy"=19 ,       "boy" =19, 
             "father"=4 ,     "mother" =0, 
             "boy"=19 ,        "boy"=19 , 
             "father"=4 ,     "mother" =0, 
             "boy"=19 ,        "boy" =19, 
             "father"=4 ,     "mother" =0, 
             "boy"=19 ,        "boy" =19, 
             "father" =4,     "mother"=0 
 
  )




myType3_color <- c("IVF-frozen"="pink" ,     "IVF-frozen"="pink" , 
             "IVF-frozen"="pink" ,     "IVF-frozen"="pink" ,   "IVF-frozen"="pink" , 
             "IVF-frozen"="pink" ,     "IVF-frozen"="pink" , 
             "IVF-frozen"="pink" ,     "IVF-frozen"="pink" , 
             "IVF-frozen" ="pink",     "IVF-frozen" ="pink", 
             "IVF-frozen" ="pink",     "IVF-frozen" ="pink", 
             "ICSI-fresh"="blue" ,     "ICSI-fresh"="blue",  "ICSI-fresh"="blue", 
             "ICSI-fresh"="blue" ,     "ICSI-fresh"="blue", 
             "ICSI-fresh"="blue" ,     "ICSI-fresh"="blue", 
             "ICSI-fresh"="blue" ,     "ICSI-fresh"="blue", 
             "ICSI-fresh" ="blue",     "ICSI-fresh"="blue",
             "ICSI-fresh"="blue" ,     "ICSI-fresh"="blue", 
             "ICSI-frozen"="skyblue" ,     "ICSI-frozen"="skyblue",
             "ICSI-frozen" ="skyblue",     "ICSI-frozen"="skyblue"  
)




## boy=19, girl=17,  father=4, mother=0
## NC=cyan,   IVF-fresh=red, IVF-frozen=pink,  ICSI-fresh=blue, ICSI-frozen=skyblue



myTreatment <- c(1:length(FileLists))  ##    This option will determine the result of unite.
myOutDir <- paste(myOutDir, "_minReads=", myMinReads, sep="")
if( ! file.exists(myOutDir) ) { dir.create(myOutDir, recursive = TRUE) }
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
library(ggbiplot)


continue_on_error <- function()
{
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set
'options(error=continue_on_error())'")
}

# This is the key option
options(error=continue_on_error) 





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











# read the files to a methylRawList object: myobj
myobj=methRead(FileLists,
               sample.id=mySampleID,
               assembly="hg38",
               treatment=myTreatment,
               context="CpG",
               pipeline = "bismarkCoverage",
               mincov = myMinReads,       ## >= n
               header = FALSE
)


sink( file=paste(myOutDir, "1-all-the-files.txt", sep="/") )
  print(FileLists)
  print("#########################")
  print("#########################")
  print(myobj)
sink()





pdf( file=paste(myOutDir, "2-MethylationStats-the-files.pdf", sep="/")  )
for( i in c(1:length(FileLists)) ) {
  getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()


 

pdf( file=paste(myOutDir, "3-CoverageStats-the-files.pdf", sep="/")  )
for( i in c(1:length(FileLists)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE )
}
dev.off()



myOutDir2 = paste(myOutDir, "/",  "sitesCov_", myMinReads,  sep="")
myOutDir3 = paste(myOutDir, "/",  "bedGraph_", myMinReads,  sep="")

if( ! file.exists(myOutDir2) ) { dir.create(myOutDir2, recursive = TRUE) }
if( ! file.exists(myOutDir3) ) { dir.create(myOutDir3, recursive = TRUE) }

for( i in c(1:length(FileLists)) ) {
  mytemp1 <- getData( myobj[[i]] )
  mytemp2 <- mytemp1[,6]/mytemp1[,5]
  mytemp3 <- cbind(mytemp1[,c(1:3)], mytemp2)
  
  write.table(myobj[[i]], file = paste(myOutDir2, "/",  "reads", myMinReads, "_",  mySampleID[[i]],  sep=""), 
              append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
  write.table(mytemp3, file = paste(myOutDir3, "/",  "reads", myMinReads, "_",  mySampleID[[i]], ".bedGraph", sep=""), 
              append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")
}



#Filtering samples based on read coverage
#It might be useful to filter samples based on coverage. Particularly, 
#if our samples are suffering from PCR bias it would be useful to discard bases with very high read coverage. 
#Furthermore, we would also like to discard bases that have low read coverage, a high enough read coverage will 
#increase the power of the statistical tests. The code below filters a methylRawList and discards bases that have 
#coverage below 10X and also discards the bases that have more than 99.9th percentile of coverage in each sample.
## ------------------------------------------------------------------------
filtered.myobj=filterByCoverage(myobj,  lo.count=myMinReads,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.9)
sink( file=paste(myOutDir, "4A-all-the-files.txt", sep="/") )
print(filtered.myobj)
sink()


discarded.filtered.myobj=filterByCoverage(myobj,  lo.count=myMinReads,  lo.perc=99.9,  hi.count=NULL, hi.perc=NULL)
sink( file=paste( myOutDir, "4B-all-the-files.discarded.txt", sep="/") )
print(discarded.filtered.myobj)
sink()

## myobj  and filtered.myobj, both of them will be used to further analysis.


sink( file=paste(myOutDir, "5A-dimensions-myobj.txt", sep="/")  )
for( i in c(1:length(FileLists)) ) {
  print(    mySampleID[[i]]  )
  print(   dim(myobj[[i]])  )
}
sink()


sink( file=paste(myOutDir, "5B-dimensions-filteredMyobj.txt", sep="/")  )
for( i in c(1:length(FileLists)) ) {
  print(    mySampleID[[i]]  )  
  print(    dim(filtered.myobj[[i]]) )  
}
sink()




######################################################################################################################################################
## myobj  and filtered.myobj, both of them will be used to further analysis.
######################################################################################################################################################
## Merging samples
## Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. 
## This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses).
meth=unite(myobj, destrand=FALSE)
head(meth)
dim(meth)


sink(file=paste(myOutDir, "6-dimensions.txt", sep="/") )
print(dim(meth))
sink()

pdf( file=paste(myOutDir, "6-getCorrelation.pdf", sep="/")  )
#getCorrelation(meth, plot=TRUE, method="pearson" )
#getCorrelation(meth, plot=TRUE, method="kendall" )
#getCorrelation(meth, plot=TRUE, method="spearman" )
dev.off()



pdf( file=paste(myOutDir, "7-clusterSamples.pdf", sep="/") , width=5, height=20  )
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="ward", plot=TRUE)
clusterSamples(meth, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



pdf( file=paste(myOutDir, "8-PCASamples-importance-of-components.pdf", sep="/")  )
PCASamples(meth, screeplot=TRUE)
dev.off()


pdf( file=paste(myOutDir, "9-PCASamples.pdf", sep="/")  )
PCASamples(meth, adj.lim=c(5, 5))
dev.off()





######################################################################################################################################################
write.table(meth, file = paste(myOutDir,"10_all_merged_overlap.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

meth2 <- getData(meth) 
meth2 <- meth2[,-c(1:4)] 
dim(meth)
dim(meth2)
meth2[1:10,]

metMatrix <- matrix(ncol=nrow(meth2), nrow=length(mySampleID) )
rownames(metMatrix) = mySampleID
for( i in c(1:length(FileLists)) ) {
  i2 = (i-1)*3+1
  i3 = (i-1)*3+2
  metMatrix[i,] = meth2[,i3]/meth2[,i2]
}

metMatrix[,1:10]
dim(metMatrix)
write.table(metMatrix , file = paste(myOutDir,"12A_me-ratio-Matrix.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



mat=percMethylation(meth)
dim(mat)

write.table(mat , file = paste(myOutDir,"12B_me-ratio-Matrix.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")


 
metMatrix[is.na(metMatrix)] <- 0

write.table(metMatrix , file = paste(myOutDir,"12C_me-ratio-Matrix_NAby0.txt",  sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")





#####################################################################################################
myPCA1 <- prcomp(metMatrix)
names(myPCA1)

sink( file = paste(myOutDir,"13A_PCA.txt",  sep="/") )
print(myPCA1)
sink()

sink( file = paste(myOutDir,"13B_PCA-summary.txt",  sep="/") )
summary(myPCA1)
sink()

sink( file = paste(myOutDir,"13C_PCA-all.txt",  sep="/") )
print("####################### myPCA1$sdev #########################")
print(myPCA1$sdev)
print("####################### myPCA1$rotation #########################")
print(myPCA1$rotation)
print("####################### myPCA1$center #########################")
print(myPCA1$center)
print("####################### myPCA1$scale #########################")
print(myPCA1$scale)
print("####################### myPCA1$x #########################")
print(myPCA1$x)
sink()


pdf( file=paste(myOutDir, "13D-PCA-info.pdf", sep="/")  )
plot(myPCA1, type="lines")
fviz_eig(myPCA1)
dev.off() 


myPCA1_matrix <- myPCA1$x
dim(myPCA1_matrix)


open3d()
plot3d(myPCA1_matrix[,1:3] , col=myType1, size=5)
rgl.postscript(paste(myOutDir, "14A-PCA-3D.pdf", sep="/"),"pdf")  



my_fviz_pca_ind1 <- fviz_pca_ind(myPCA1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
my_fviz_pca_ind2 <- fviz_pca_ind(myPCA1,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE
)
my_fviz_pca_ind3 <- fviz_pca_ind(myPCA1,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 addEllipses = TRUE, # Concentration ellipses
                                 ellipse.type = "confidence",
                                 legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)
my_fviz_pca_ind4 <- fviz_pca_ind(myPCA1,
                                 col.ind =  as.factor(myType3), # color by groups
                                 #palette = c("#00AFBB",  "#FC4E07"),
                                 #legend.title = "Groups",
                                 repel = TRUE, 
                                 label = "none", 
                                 alpha.ind = 1
)

pdf( file=paste(myOutDir, "14B-PCA-2D-1.pdf", sep="/")  )
print(my_fviz_pca_ind1)
dev.off() 
svg(file=paste(myOutDir, "14B-PCA-2D-1.svg", sep="/") )
print(my_fviz_pca_ind1)
dev.off() 

pdf( file=paste(myOutDir, "14B-PCA-2D-2.pdf", sep="/")  )
print(my_fviz_pca_ind2)
dev.off() 
svg(file=paste(myOutDir, "14B-PCA-2D-2.svg", sep="/") )
print(my_fviz_pca_ind2)
dev.off() 

pdf( file=paste(myOutDir, "14B-PCA-2D-3.pdf", sep="/")  )
print(my_fviz_pca_ind3)
dev.off() 
svg(file=paste(myOutDir, "14B-PCA-2D-3.svg", sep="/") )
print(my_fviz_pca_ind3)
dev.off() 

pdf( file=paste(myOutDir, "14B-PCA-2D-4.pdf", sep="/")  )
print(my_fviz_pca_ind4)
dev.off() 
svg(file=paste(myOutDir, "14B-PCA-2D-4.svg", sep="/") )
print(my_fviz_pca_ind4)
dev.off() 
 
  


dataframeA  <- data.frame( as.data.frame(myPCA1_matrix), myType2, myType3) 
dataframeA 

 
FigureTemp1 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
               geom_point(size=3, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
               scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
               MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
               guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5)))  
MySaveGgplot2_1(ggplot2Figure1=FigureTemp1,  path1=myOutDir, fileName1="15A-PCA-PC1-PC2",  height1=4,  width1=6)
 


FigureTemp2 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=3, alpha=0.5  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp2,  path1=myOutDir, fileName1="15B-PCA-PC1-PC2-alpha",  height1=4,  width1=6)


FigureTemp3 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=1, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp3,  path1=myOutDir, fileName1="15C-PCA-PC1-PC2-smallDot",  height1=4,  width1=6)



FigureTemp4 <- ggplot( data = dataframeA, aes(x = PC1, y = PC2, shape=as.factor(myType2), color=as.factor(myType3) )) + 
  geom_point(size=5, alpha=1  ) + xlab("PC1") +   ylab("PC2") +   
  scale_colour_manual(values=myType3_color) + scale_shape_manual(values=myType2_shape) + 
  MyTheme_1( textSize=14,  hjust1=NULL, vjust1=NULL,  angle1=NULL) +
  guides( colour = guide_legend(override.aes = list(size=5)) ,  shape = guide_legend(override.aes = list(size=5))) 
MySaveGgplot2_1(ggplot2Figure1=FigureTemp4,  path1=myOutDir, fileName1="15D-PCA-PC1-PC2-big",  height1=4,  width1=6)




#####################################################################################
myCluster1 <- clara(metMatrix, 2)
names(myCluster1)
sink( file = paste(myOutDir,"16A_myCluster1-clara-2classes.txt",  sep="/") )
myCluster1$clusinfo
myCluster1$clustering
myCluster1$sample
myCluster1$i.med
myCluster1$objective  
myCluster1$diss  
myCluster1$call  
myCluster1$silinfo  
sink() 
   



myCluster1 <- clara(metMatrix, 3)
names(myCluster1)
sink( file = paste(myOutDir,"16A_myCluster1-clara-3classes.txt",  sep="/") )
myCluster1$clusinfo
myCluster1$clustering
myCluster1$sample
myCluster1$i.med
myCluster1$objective  
myCluster1$diss  
myCluster1$call  
myCluster1$silinfo  
sink() 



myCluster1 <- clara(metMatrix, 4)
names(myCluster1)
sink( file = paste(myOutDir,"16A_myCluster1-clara-4classes.txt",  sep="/") )
myCluster1$clusinfo
myCluster1$clustering
myCluster1$sample
myCluster1$i.med
myCluster1$objective  
myCluster1$diss  
myCluster1$call  
myCluster1$silinfo  
sink() 



myCluster1 <- clara(metMatrix, 5)
names(myCluster1)
sink( file = paste(myOutDir,"16A_myCluster1-clara-5classes.txt",  sep="/") )
myCluster1$clusinfo
myCluster1$clustering
myCluster1$sample
myCluster1$i.med
myCluster1$objective  
myCluster1$diss  
myCluster1$call  
myCluster1$silinfo  
sink() 



myCluster2 <- fanny(metMatrix, 2)
names(myCluster2)
sink( file = paste(myOutDir,"16B_myCluster2-fanny-2classes.txt",  sep="/") )
myCluster2$membership
myCluster2$coeff
myCluster2$memb.exp
myCluster2$clustering
myCluster2$k.crisp
myCluster2$convergence
myCluster2$objective  
myCluster2$diss  
myCluster2$silinfo  
sink() 


 
 



