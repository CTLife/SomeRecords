##############################################################################################################################################################################################
## Flanks_Dis_DRACH                     
## Author: Yong Peng, yongp@outlook.com 
## Version 0.0.1, March 1st, 2020.
## Run "Rscript Flanks_Dis_DRACH.R -h" to get some help.
##############################################################################################################################################################################################





##############################################################################################################################################################################################
suppressPackageStartupMessages( library(optparse) )  ## To run the script in command lines.

getParameters_f <- function() {
  option_list_Local <- list(   ## Options list with associated default value.  
      optparse::make_option(opt_str=c("-A", "--inputA"),
      default="1-rawFiles/WT1_vs_Input/minus.T.txt",
      type="character",   dest="inputA",
      help="Input file with all information. [default: %default]."),

      optparse::make_option(opt_str=c("-B", "--inputB"),
      default="12-lineNumbers-5%/WT1_vs_Input/DRACH.minus.T/withRIP_FP_DRACH.txt",
      type="character",   dest="inputB",
      help="Input file for FP sites. [default: %default]."),

      optparse::make_option(opt_str=c("-C", "--inputC"),
      default="12-lineNumbers-5%/WT1_vs_Input/DRACH.minus.T/withRIP_TP_DRACH.txt",
      type="character",   dest="inputC",
      help="Input file for TP sites. [default: %default]."),

      optparse::make_option(opt_str=c("-O", "--outDir"),
      default="22-lineNumbers-5percent/WT1_vs_Input/DRACH.minus.T.withRIP",
      type="character",   dest="outDir",
      help="Path to the directory containing all the analysis results. [default: %default]."),

      optparse::make_option(opt_str=c("-f", "--flanks"),
      default=100,
      type="integer",   dest="flanks",
      help="Number of upstream (downstream) as background. [default: %default].")          
)
  
  ## Now parse the command line to check which option is given and get associated values.
  parser_Local <- optparse::OptionParser(usage="usage: Rscript %prog [options]",
      option_list=option_list_Local, 
      description="Flanks_Dis_DRACH.R, Version 0.0.1, March 1st, 2020.",                             
      epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>."
  )
  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)
}
##############################################################################################################################################################################################


 


##############################################################################################################################################################################################
opt_g = getParameters_f()  

inputA_g      <- opt_g$inputA
inputB_g      <- opt_g$inputB
inputC_g      <- opt_g$inputC
outDir_g      <- opt_g$outDir
flanks_g      <- as.integer(opt_g$flanks)
 
rm(getParameters_f)  
rm(opt_g)  

options(digits=10)
continue_on_error_g <- function() {
    print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'. " )
}
options( error=continue_on_error_g )  ## This option is very important.

print("####################")
print(inputA_g)
print(inputB_g)
print(inputC_g)
print(outDir_g)
print(flanks_g)
print("####################")
##############################################################################################################################################################################################





############################################################################################################################################################################################## Load all required packages
library(tidyverse)
library(stringr)
library(ggplot2)
library(gsubfn)
##############################################################################################################################################################################################






############################################################################################################################################################################################## read input data
# inputA  = "1-rawFiles/WT1_vs_Input/minus.T.txt"
# inputFP = "12-lineNumbers-5%/WT1_vs_Input/DRACH.minus.T/withRIP_FP_DRACH.txt"
# inputTP = "12-lineNumbers-5%/WT1_vs_Input/DRACH.minus.T/withRIP_TP_DRACH.txt"
# outDir  = "22-lineNumbers-5%/WT1_vs_Input/DRACH.minus.T.withRIP"

dir.create(path = outDir_g, showWarnings = TRUE, recursive = TRUE, mode = "0777")

dataFrameA  = read.table(file = inputA_g, header = TRUE, sep = "\t", quote = "\"'", dec = ".", numerals = "no.loss",  na.strings = "NA", colClasses = NA,              
                         nrows = -1, skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE )
dataFrameFP = read.table(file = inputB_g, header = FALSE, sep = "\t", quote = "\"'", dec = ".", numerals = "no.loss",  na.strings = "NA", colClasses = NA,              
                         nrows = -1, skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE )
dataFrameTP = read.table(file = inputC_g, header = FALSE, sep = "\t", quote = "\"'", dec = ".", numerals = "no.loss",  na.strings = "NA", colClasses = NA,              
                         nrows = -1, skip = 0, check.names = TRUE, strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE )

dim(dataFrameA)
dim(dataFrameFP)
dim(dataFrameTP)

head(dataFrameA)
head(dataFrameFP)
head(dataFrameTP)

vectorFP = as.character(dataFrameFP[,1]) 
vectorTP = as.character(dataFrameTP[,1]) 
length(vectorFP)
length(vectorTP)

vectorFP_1 = gsubfn::strapplyc(X=vectorFP, pattern="^(\\d+):chr" ) 
vectorFP_1 = as.numeric( as.vector(unlist(vectorFP_1)) )

vectorTP_1 = gsubfn::strapplyc(X=vectorTP, pattern="^(\\d+):chr" ) 
vectorTP_1 = as.numeric( as.vector(unlist(vectorTP_1)) )
##############################################################################################################################################################################################





############################################################################################################################################################################################## generate background matrix 
matrixFP = matrix(nrow=length(vectorFP_1) , ncol=flanks_g*2)
matrixTP = matrix(nrow=length(vectorTP_1) , ncol=flanks_g*2)

dim(matrixFP)
dim(matrixTP)

j=1
for(i in vectorFP_1 ) {
  upstream = c((i-flanks_g-1):(i-2))
  downstream = c((i+0):(i+flanks_g-1))
  index1 = c(upstream, downstream)
  index1[index1<1] = 1
  index1[index1>length(dataFrameA$tumor_var_freq)] = length(dataFrameA$tumor_var_freq)
  matrixFP[j,] = as.character( dataFrameA$tumor_var_freq[index1] )
  j = j+1
}

j=1
for(i in vectorTP_1 ) {
  upstream = c((i-flanks_g-1):(i-2))
  downstream = c((i+0):(i+flanks_g-1))
  index1 = c(upstream, downstream)
  index1[index1<1] = 1
  index1[index1>length(dataFrameA$tumor_var_freq)] = length(dataFrameA$tumor_var_freq)
  matrixTP[j,] = as.character( dataFrameA$tumor_var_freq[index1] )
  j = j+1
}


dim(matrixFP)
dim(matrixTP)

write.table(x=matrixFP, file = paste(outDir_g, "1-matrixFP.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")
write.table(x=matrixTP, file = paste(outDir_g, "1-matrixTP.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE, qmethod = c("escape", "double"), fileEncoding = "")

matrixFP_1 = str_replace(string=matrixFP, pattern="%", replacement="")
matrixTP_1 = str_replace(string=matrixTP, pattern="%", replacement="")
matrixFP_1 = as.numeric(matrixFP_1 )
matrixTP_1 = as.numeric(matrixTP_1 )
matrixFP_1 = matrix(data=matrixFP_1 , nrow=length(vectorFP_1) , ncol=flanks_g*2)
matrixTP_1 = matrix(data=matrixTP_1 , nrow=length(vectorTP_1) , ncol=flanks_g*2)

min_row = min( nrow(matrixFP), nrow(matrixTP) )
matrixFP_1_random = matrix(nrow=min_row, ncol=flanks_g*2)
matrixTP_1_random = matrix(nrow=min_row, ncol=flanks_g*2)

if(nrow(matrixFP) >= nrow(matrixTP)) { 
    index_random = sample(x=c(1:nrow(matrixFP)), size=nrow(matrixTP) )  
    matrixFP_1_random = matrixFP_1[index_random,]
    matrixTP_1_random = matrixTP_1
}
if(nrow(matrixFP) <  nrow(matrixTP)) { 
    index_random = sample(x=c(1:nrow(matrixTP)), size=nrow(matrixFP) )  
    matrixFP_1_random = matrixFP_1
    matrixTP_1_random = matrixTP_1[index_random,]
}



dim(matrixFP_1)
dim(matrixTP_1)
dim(matrixFP_1_random)
dim(matrixTP_1_random)

write.table(x=matrixFP_1, file = paste(outDir_g, "2a-matrixFP_1.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE, fileEncoding = "")
write.table(x=matrixTP_1, file = paste(outDir_g, "2b-matrixTP_1.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE, fileEncoding = "")
write.table(x=matrixFP_1_random, file = paste(outDir_g, "2c-matrixFP_1_random.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE, fileEncoding = "")
write.table(x=matrixTP_1_random, file = paste(outDir_g, "2d-matrixTP_1_random.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE, fileEncoding = "")
##############################################################################################################################################################################################





############################################################################################################################################################################################## dataframe as input of ggplot2 
xAxis1  = c(-flanks_g:-1, 1:flanks_g)
xAxisFP = rep(xAxis1, nrow(matrixFP_1))
xAxisTP = rep(xAxis1, nrow(matrixTP_1))
xAxisFP_random = rep(xAxis1, nrow(matrixFP_1_random))
xAxisTP_random = rep(xAxis1, nrow(matrixTP_1_random))

yAxisFP = as.vector(t(matrixFP_1))
yAxisTP = as.vector(t(matrixTP_1))
yAxisFP_random = as.vector(t(matrixFP_1_random))
yAxisTP_random = as.vector(t(matrixTP_1_random))

typeFP  = rep("FP",length(as.vector(t(matrixFP_1))))
typeTP  = rep("TP",length(as.vector(t(matrixTP_1))))
typeFP_random  = rep("FP",length(as.vector(t(matrixFP_1_random))))
typeTP_random  = rep("TP",length(as.vector(t(matrixTP_1_random))))

xAxis2 = c(xAxisFP  ,  xAxisTP)
yAxis2 = c(yAxisFP ,   yAxisTP)
type2  = c(typeFP  ,  typeTP)
 
xAxis2_random = c(xAxisFP_random  ,  xAxisTP_random)
yAxis2_random = c(yAxisFP_random ,   yAxisTP_random)
type2_random  = c(typeFP_random  ,   typeTP_random)

print("####################")
length(xAxis2)
length(yAxis2)
length(type2)
print("####################")
length(xAxis2_random)
length(yAxis2_random)
length(type2_random)
print("####################")

yAxis2[yAxis2>10] = 10 
yAxis2_random[yAxis2_random>10] = 10 

dataframeA         <- data.frame( xAxis=xAxis2 , yAxis=yAxis2 , sampleType=type2 )
dataframeA_random  <- data.frame( xAxis=xAxis2_random, yAxis=yAxis2_random, sampleType=type2_random )
##############################################################################################################################################################################################





############################################################################################################################################################################################## for figures
MyTheme_g <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  ggplot2::theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),     
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    
    ## aspect.ratio = 1,       
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),        
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),        
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),        
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),                                                          
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),        
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),      
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm. 
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## lines along axes (element_line; inherits from line). 
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                ## extra space added around legend (unit). 
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	## background underneath legend keys. 
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	            ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")    
    legend.justification = "center",      	    ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)   
    legend.box           = NULL, 	            ## arrangement of multiple legends ("horizontal" or "vertical")  
    legend.box.just      = NULL, 	            ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")       
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),     ## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)                                     
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                            ## margin around facet panels (unit)   
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                            ## grid lines (element_line; inherits from line)   
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,       ## minor grid lines (element_line; inherits from panel.grid)   
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,       ## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,       ## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,       ## horizontal minor grid lines (element_line; inherits from panel.grid.minor)    
    plot.background  = element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                            ## background of the entire plot (element_rect; inherits from rect)   
    plot.title       = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)   
    plot.margin      = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)   
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	      ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	      ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	      ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 


MySaveGgplot2_g <- function(path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  #EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  #if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  ggplot2::ggsave(filename=paste(SVG1, "/", fileName1, ".svg", sep=""),  plot = last_plot(), device = "svg",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 1000, limitsize = FALSE)
  ggplot2::ggsave(filename=paste(PDF1, "/", fileName1, ".pdf", sep=""),  plot = last_plot(), device = "pdf",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 1000, limitsize = FALSE)
  #ggplot2::ggsave(filename=paste(EPS1, "/", fileName1, ".eps", sep=""),  plot = last_plot(), device =cairo_ps, path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 1000, limitsize = FALSE)
  ggplot2::ggsave(filename=paste(PNG1, "/", fileName1, ".png", sep=""),  plot = last_plot(), device = "png",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 1000, limitsize = FALSE )
}
##############################################################################################################################################################################################





############################################################################################################################################################################################## 
## scatter plot
ggplot(data=dataframeA_random, aes(x=xAxis, y=yAxis, fill=sampleType, colour=sampleType) )   +  
    xlab("Position") + ylab("Mutation Frequency (%)") +   geom_point(alpha=0.5, size=0.5 ) + 
    MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
MySaveGgplot2_g(path1=outDir_g, fileName1="1a_random_all",  height1=6,  width1=10)

dataframeA_FP_random  <- data.frame( xAxis=xAxisFP_random, yAxis=yAxisFP_random, sampleType=typeFP_random )
dataframeA_FP_random$yAxis[dataframeA_FP_random$yAxis>10] = 10
ggplot(data=dataframeA_FP_random, aes(x=xAxis, y=yAxis, fill=sampleType, colour=sampleType) )   +  
  xlab("Position") + ylab("Mutation Frequency (%)") +   geom_point(alpha=0.8, size=0.5, colour="red") + 
  MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
MySaveGgplot2_g(path1=outDir_g, fileName1="1b_random_FP",  height1=6,  width1=10)

dataframeTP_random  <- data.frame( xAxis=xAxisTP_random, yAxis=yAxisTP_random, sampleType=typeTP_random )
dataframeTP_random$yAxis[dataframeTP_random$yAxis>10] = 10
ggplot(data=dataframeTP_random, aes(x=xAxis, y=yAxis, fill=sampleType, colour=sampleType) )   +  
  xlab("Position") + ylab("Mutation Frequency (%)") +   geom_point(alpha=0.5, size=0.5, colour="red") + 
  MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
MySaveGgplot2_g(path1=outDir_g, fileName1="1c_random_TP",  height1=6,  width1=10)
## coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE)



## scatter plot
ggplot(data=dataframeA_random, aes(x=xAxis, y=yAxis, fill=sampleType, colour=sampleType) )   +  
    xlab("Position") + ylab("Mutation Frequency (%)") +   geom_point(alpha=0.5, size=0.5 ) + 
    MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) + ylim(0,5) 
MySaveGgplot2_g(path1=outDir_g, fileName1="2a_random_all",  height1=6,  width1=10)

dataframeA_FP_random  <- data.frame( xAxis=xAxisFP_random, yAxis=yAxisFP_random, sampleType=typeFP_random )
dataframeA_FP_random$yAxis[dataframeA_FP_random$yAxis>10] = 10
ggplot(data=dataframeA_FP_random, aes(x=xAxis, y=yAxis, fill=sampleType, colour=sampleType) )   +  
  xlab("Position") + ylab("Mutation Frequency (%)") +   geom_point(alpha=0.8, size=0.5, colour="red") + 
  MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  + ylim(0,5)   
MySaveGgplot2_g(path1=outDir_g, fileName1="2b_random_FP",  height1=6,  width1=10)

dataframeTP_random  <- data.frame( xAxis=xAxisTP_random, yAxis=yAxisTP_random, sampleType=typeTP_random )
dataframeTP_random$yAxis[dataframeTP_random$yAxis>10] = 10
ggplot(data=dataframeTP_random, aes(x=xAxis, y=yAxis, fill=sampleType, colour=sampleType) )   +  
  xlab("Position") + ylab("Mutation Frequency (%)") +   geom_point(alpha=0.5, size=0.5, colour="red") + 
  MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  + ylim(0,5)   
MySaveGgplot2_g(path1=outDir_g, fileName1="2c_random_TP",  height1=6,  width1=10)
## coord_cartesian(xlim = NULL, ylim = NULL, expand = TRUE)
##############################################################################################################################################################################################





##############################################################################################################################################################################################  
## y-axis is number  
ggplot(dataframeA_random, aes(x=yAxis, colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
      MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
MySaveGgplot2_g(path1=outDir_g, fileName1="4a_random_number",  height1=6,  width1=10)

ggplot(dataframeA_random, aes(x=yAxis, colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
       MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   + 
       coord_cartesian(xlim = c(0.1, 5),  ylim = c(0, 2000), expand = FALSE)  
MySaveGgplot2_g(path1=outDir_g, fileName1="4b_random_number",  height1=6,  width1=10)

ggplot(dataframeA_random, aes(x=yAxis, colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
       MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    + 
       coord_cartesian(xlim = c(0.1, 5),  ylim = c(0, 1000), expand = FALSE)    
MySaveGgplot2_g(path1=outDir_g, fileName1="4c_random_number",  height1=6,  width1=10)

ggplot(dataframeA_random, aes(x=yAxis, colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
       MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   + 
       coord_cartesian(xlim = c(0.1, 3),  ylim = c(0, 2000), expand = FALSE)  
MySaveGgplot2_g(path1=outDir_g, fileName1="4d_random_number",  height1=6,  width1=10)

ggplot(dataframeA_random, aes(x=yAxis, colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
       MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )    + 
       coord_cartesian(xlim = c(0.1, 3),  ylim = c(0, 1000), expand = FALSE)    
MySaveGgplot2_g(path1=outDir_g, fileName1="4e_random_number",  height1=6,  width1=10)
##############################################################################################################################################################################################





##############################################################################################################################################################################################  
## y-axis is density 
ggplot(dataframeA_random, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
MySaveGgplot2_g(path1=outDir_g, fileName1="5a_random_density",  height1=6,  width1=10)

ggplot(dataframeA_random, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) +
   coord_cartesian(xlim = c(0.1, 5),  ylim = c(0, 0.8), expand = FALSE)  
MySaveGgplot2_g(path1=outDir_g, fileName1="5b_random_density",  height1=6,  width1=10)

ggplot(dataframeA_random, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  +
   coord_cartesian(xlim = c(0.1, 5),  ylim = c(0, 0.5), expand = FALSE) 
MySaveGgplot2_g(path1=outDir_g, fileName1="5c_random_density",  height1=6,  width1=10)




ggplot(dataframeA, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )   
MySaveGgplot2_g(path1=outDir_g, fileName1="6a_density",  height1=6,  width1=10)

ggplot(dataframeA, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) +
   coord_cartesian(xlim = c(0.1, 5),  ylim = c(0, 0.6), expand = FALSE)  
MySaveGgplot2_g(path1=outDir_g, fileName1="6b_density",  height1=6,  width1=10)

ggplot(dataframeA, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  +
   coord_cartesian(xlim = c(0.1, 5),  ylim = c(0, 0.3), expand = FALSE) 
MySaveGgplot2_g(path1=outDir_g, fileName1="6c_density",  height1=6,  width1=10)

ggplot(dataframeA, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 ) +
   coord_cartesian(xlim = c(0.1, 2.5),  ylim = c(0, 0.6), expand = FALSE)  
MySaveGgplot2_g(path1=outDir_g, fileName1="6d_density",  height1=6,  width1=10)

ggplot(dataframeA, aes(x=yAxis, stat(density), colour = sampleType)) + geom_freqpoly(binwidth = 0.1) + 
   MyTheme_g( hjust1=NULL, vjust1=NULL,  angle1=NULL,   textSize=14 )  +
   coord_cartesian(xlim = c(0.1, 2.5),  ylim = c(0, 0.3), expand = FALSE) 
MySaveGgplot2_g(path1=outDir_g, fileName1="6e_density",  height1=6,  width1=10)
##############################################################################################################################################################################################










