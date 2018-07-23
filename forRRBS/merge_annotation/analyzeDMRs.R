#! /usr/bin/env Rscript
###########################################################################################################
### analyzeDMRs: Analyze raw DMRs from methylKit, such as merge, format and annotation.
### Author: Yong Peng
### Run "./analyzeDMRs.R  -h"  or "Rscript analyzeDMRs.R  -h" to get some help.
###########################################################################################################





## Part 1: To run the script in command lines.
#################################################
suppressPackageStartupMessages( library(optparse) )      

getParameters_f <- function() {
  option_list_Local <- list(   # Options list with associated default value.   
      optparse::make_option(opt_str=c("-f", "--file"),
			default="5_DMR/2E_AlldiffMesites_q0.05_diff0.txt",
			type="character",   dest="file",
			help="Path to the file with DMRs or DMCs information from methylKit. [default: %default]."),

      optparse::make_option(opt_str=c("-o", "--outDir"),
			default="5_DMR_Annotation",
			type="character",   dest="outDir",
			help="Path or dir name for output files. [default: %default]."),

      optparse::make_option(opt_str=c("-rg", "--RefGenome"),
			default="hg38",
			type="character",   dest="RefGenome",
			help="One of 'hg38' and 'mm10'.   [default: %default].") , 

      optparse::make_option(opt_str=c("-name", "--RegionName"),
			default="NA",
			type="character",   dest="RegionName",
			help="The prefix of name (the 4th column in bed file) of each genomic region.   [default: %default].") 
  )

  # now parse the command line to check which option is given and get associated values
  parser_Local <- optparse::OptionParser(usage="usage: %prog [options]",
		option_list=option_list_Local, 
		description="analyzeDMRs: Analyze raw DMRs from methylKit, such as merge, format and annotation. July 30th, 2018.",                             
		epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>"
  )

  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)   
}


opt_g = getParameters_f()
inFile_g     <- opt_g$file
outDir_g     <- opt_g$outDir
RefGenome_g  <- opt_g$RefGenome
RegionName_g <- opt_g$RegionName

options(digits=10)
continue_on_error <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())' " )
}
options( error=continue_on_error )  ## This option is very important.
#############################################





## Part 2: load all the required packages.
#############################################
suppressPackageStartupMessages( library(GenomicFeatures) ) 
suppressPackageStartupMessages( library(GenomicRanges) ) 
suppressPackageStartupMessages( library(clusterProfiler) ) 
suppressPackageStartupMessages( library(ReactomePA) ) 
suppressPackageStartupMessages( library(ChIPseeker) ) 
suppressPackageStartupMessages( library(DOSE) ) 
suppressPackageStartupMessages( library(ggplot2) ) 
suppressPackageStartupMessages( library(topGO) ) 
suppressPackageStartupMessages( library(KEGG.db) ) 
my_txdb_g  <- ""
my_orgdb_g <- ""
if(RefGenome_g == "hg38") {
    suppressPackageStartupMessages( library(TxDb.Hsapiens.UCSC.hg38.knownGene) ) 
    suppressPackageStartupMessages( library(org.Hs.eg.db) )  
    my_txdb_g  <- TxDb.Hsapiens.UCSC.hg38.knownGene
    my_orgdb_g <- org.Hs.eg.db
}
if(RefGenome_g == "mm10") {
    suppressPackageStartupMessages( library(TxDb.Mmusculus.UCSC.mm10.knownGene) ) 
    suppressPackageStartupMessages( library(org.Mm.eg.db) )  
    my_txdb_g  <- TxDb.Mmusculus.UCSC.mm10.knownGene
    my_orgdb_g <- org.Mm.eg.db
}
#############################################





## Part 3: Select the DMRs or DMCs.
#############################################
mySelectDiffMe <- function( path_temp3, file_temp3, name3 ) {
    myOutPath = paste(path_temp3, "1-SelectedDiffMe", sep="/")
    if( ! file.exists(myOutPath) )  { dir.create(myOutPath, recursive = TRUE) }
    nameOfRegion = name3
    sink( file=paste(myOutPath, "runLog.txt",  sep="/" ) )
    matrix1 = read.table(file= file_temp3, header = TRUE, sep = "\t", quote = "\"'", dec = ".")
    dim(matrix1)
    colnames(matrix1)
    myQvalue = matrix1$qvalue
    myDiff   = matrix1$meth.diff   
    print("#####################################################")
    print("## Number of DMCs for six conditions: ")
    print(nrow( matrix1[( (myQvalue<0.001) & (abs(myDiff)>10) ),]  ) )
    print(nrow( matrix1[( (myQvalue<0.01)  & (abs(myDiff)>10) ),]  ) )
    print(nrow( matrix1[( (myQvalue<0.05)  & (abs(myDiff)>10) ),]  ) )
    print(nrow( matrix1[( (myQvalue<0.001) & (abs(myDiff)>5) ),]  ) )
    print(nrow( matrix1[( (myQvalue<0.01)  & (abs(myDiff)>5) ),]  ) )
    print(nrow( matrix1[( (myQvalue<0.05)  & (abs(myDiff)>5) ),]  ) )
    print("#####################################################")

    ## Select the DMCs.
    diff1_all = matrix1[( (myQvalue<0.05)  & (abs(myDiff)>5) ),]
    print("#####################################################")
    print("## Dimension of diff1_all : ")
    print(dim(diff1_all)  )
    print("#####################################################")
    diff1_name = c()
    for(i in c(1:nrow(diff1_all)) ) {
      diff1_name[i] = paste(nameOfRegion, i, sep="_"  )     
    }

    ############### for ChIPseeker
    diff1_ChIPseeker = cbind(diff1_all[,1:4],  diff1_name, rep(1, nrow(diff1_all)),   diff1_all[,5:7] )   
    colnames(diff1_ChIPseeker) = c("chr", "start",   "end",  "strand", "name",  "strength",  "pvalue",    "qvalue",  "meth.diff")
    diff1_ChIPseeker[,3] = diff1_ChIPseeker[,3] + 1
    diff1_ChIPseeker_hyper = diff1_ChIPseeker[ diff1_ChIPseeker$meth.diff > 0, ]
    diff1_ChIPseeker_hypo  = diff1_ChIPseeker[ diff1_ChIPseeker$meth.diff < 0, ]
    print("#####################################################")
    print("diff1_ChIPseeker:")
    print(dim(diff1_ChIPseeker)) 
    print(dim(diff1_ChIPseeker_hyper) )
    print(dim(diff1_ChIPseeker_hypo)) 
    print("#####################################################")

    ############### for BED
    diff1_BED = cbind(diff1_all[,1:3],  diff1_name )   
    colnames(diff1_BED) = c("chr", "start",   "end",  "name" )
    diff1_BED[,3] = diff1_BED[,3] + 1
    diff1_BED_hyper = diff1_BED[ diff1_ChIPseeker$meth.diff > 0, ]
    diff1_BED_hypo  = diff1_BED[ diff1_ChIPseeker$meth.diff < 0, ]
    print("#####################################################")
    print("diff1_BED:")
    print(dim(diff1_BED) )
    print(dim(diff1_BED_hyper) )
    print(dim(diff1_BED_hypo) )
    print("#####################################################")

    ############### for BED (not plus 1 for end site)
    diff1_BED_keptRaw = cbind(diff1_all[,1:3],  diff1_name )   
    colnames(diff1_BED_keptRaw) = c("chr", "start",   "end",  "name" )
    diff1_BED_keptRaw_hyper = diff1_BED_keptRaw[ diff1_ChIPseeker$meth.diff > 0, ]
    diff1_BED_keptRaw_hypo  = diff1_BED_keptRaw[ diff1_ChIPseeker$meth.diff < 0, ]
    print("#####################################################")
    print("diff1_BED_keptRaw:")
    print(dim(diff1_BED_keptRaw) ) 
    print(dim(diff1_BED_keptRaw_hyper))  
    print(dim(diff1_BED_keptRaw_hypo) ) 
    print("#####################################################")


    Part1_local = paste(myOutPath, "ChIPseeker",  sep="/" )
    if( ! file.exists(Part1_local) )  { dir.create(Part1_local, recursive = TRUE) }
    write.table(x=diff1_ChIPseeker,       file = paste(Part1_local, "diff1_all.txt", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
    write.table(x=diff1_ChIPseeker_hyper, file = paste(Part1_local, "diff1_hyper.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )
    write.table(x=diff1_ChIPseeker_hypo,  file = paste(Part1_local, "diff1_hypo.txt", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = TRUE  )

    Part1_local = paste(myOutPath, "BED",  sep="/" )
    if( ! file.exists(Part1_local) )  { dir.create(Part1_local, recursive = TRUE) }
    write.table(x=diff1_BED,       file = paste(Part1_local, "diff1_all.bed", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
    write.table(x=diff1_BED_hyper, file = paste(Part1_local, "diff1_hyper.bed", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
    write.table(x=diff1_BED_hypo,  file = paste(Part1_local, "diff1_hypo.bed", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )

    Part1_local =paste(myOutPath,  "BED_keptRaw",  sep="/" )
    if( ! file.exists(Part1_local) )  { dir.create(Part1_local, recursive = TRUE) }
    write.table(x=diff1_BED_keptRaw,       file = paste(Part1_local, "diff1_all.raw.bed", sep="/"),   append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
    write.table(x=diff1_BED_keptRaw_hyper, file = paste(Part1_local, "diff1_hyper.raw.bed", sep="/"), append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
    write.table(x=diff1_BED_keptRaw_hypo,  file = paste(Part1_local, "diff1_hypo.raw.bed", sep="/"),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE,  col.names = FALSE  )
    sink() 
}


if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }
mySelectDiffMe( path_temp3=outDir_g, file_temp3=inFile_g, name3=RegionName_g ) 

sink( paste(outDir_g, "1-SelectedDiffMe",  "Parameters.txt", sep="/")   )
    print( paste("inFile_g: "      , inFile_g,     sep=" ") )
    print( paste("outDir_g: "      , outDir_g,     sep=" ") )
    print( paste("RefGenome_g: "   , RefGenome_g,  sep=" ") )
    print( paste("RegionName_g: "  , RegionName_g, sep=" ") ) 
    cat("\n\n\n\n\n") 
    print( "my_txdb_g:"  )
    print(  my_txdb_g  )
    cat("\n\n\n\n\n") 
    print( "my_orgdb_g: ")  
    print(  my_orgdb_g )  
sink()   
#############################################





## Part 4: Merge the overlapped DMCs or DMRs. 
#############################################
myOutPath_1 = paste(outDir_g, "1-SelectedDiffMe", "BED", sep="/")
myOutPath_2 = paste(outDir_g, "2-Merged",   sep="/")
if( ! file.exists(myOutPath_2)   ) { dir.create(myOutPath_2, recursive = TRUE) }
myCommands_1 =  paste(  " mergePeaks  -d given ",     myOutPath_1, "/diff1_hyper.bed ",    " > ", myOutPath_2, "/diff1_hyper.txt 2>&1 ", sep="" )
myCommands_2 =  paste(  " mergePeaks  -d given ",     myOutPath_1, "/diff1_hypo.bed ",     " > ", myOutPath_2, "/diff1_hypo.txt  2>&1 ", sep="" )
system( myCommands_1 )
system( myCommands_2 )
#############################################





## Part 5: Merge the overlapped DMCs or DMRs. 
#############################################
myOutPath_3  = paste(outDir_g, "3-Formatted",   sep="/")  
if( ! file.exists(myOutPath_3)   ) { dir.create(myOutPath_3, recursive = TRUE) }
myCommands_3 =  paste(  "   perl  format.pl  -in  ", myOutPath_2,  "  -out  " , myOutPath_3,  "  -name  ", RegionName_g ) 
system( myCommands_3 )
#############################################





## Part 6: The distribution of length of genomic regions. 
#############################################

MyTheme_1_g <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
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
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),    ## line along x axis (element_line; inherits from axis.line)
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


MySaveGgplot2_1_g <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
    SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
    PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
    PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
    EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
    if( ! file.exists(SVG1) ) { dir.create(SVG1) }
    if( ! file.exists(PNG1) ) { dir.create(PNG1) }
    if( ! file.exists(PDF1) ) { dir.create(PDF1) }
    if( ! file.exists(EPS1) ) { dir.create(EPS1) }
    ggsave(filename=paste(SVG1, "/", fileName1, ".svg", sep=""),  plot = last_plot(), device = "svg",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)
    ggsave(filename=paste(PDF1, "/", fileName1, ".pdf", sep=""),  plot = last_plot(), device = "pdf",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)
    ggsave(filename=paste(EPS1, "/", fileName1, ".eps", sep=""),  plot = last_plot(), device =cairo_ps, path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)
    if( (height1<50) & (width1<50) ) {
    ggsave(filename=paste(PNG1, "/", fileName1, ".png", sep=""),  plot = last_plot(), device = "png",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)    
    }
}


## df contains two columns, the first column (cond_col=1) is sample type, the second column (val_col=2) is value. (must be).
whisk_1_g <- function(df, cond_col=1, val_col=2) {  
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



## 1 feature as type
MyBoxViolinPlot_1_f <- function(vector2,   sampleType2,  colours2,   path2,   fileName2,  title2,  xLab2,  yLab2,    height2=4,   width2=4,   Ymin2=0, Ymax2=3) { 
  vector2[vector2>Ymax2] <- Ymax2
  vector2[vector2<Ymin2] <- Ymin2
  DataFrame_Local  <- data.frame(   sampleType=sampleType2,   yAxis=vector2    ) 
  
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis ) ) +  
    geom_boxplot( outlier.shape=NA, outlier.size=0, notch=FALSE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot_NoFacet",   sep="",  collapse=NULL),  height1=height2, width1=width2-1)
  
  FigureTemp1a <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleType ) ) +  
    geom_boxplot( outlier.shape=NA, outlier.size=0, notch=FALSE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1a,  path1=path2, fileName1=paste(fileName2, "_boxPlot_NoFacet2",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp1b <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleType ) ) +  
    geom_boxplot( outlier.shape=NA, outlier.size=0, notch=FALSE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1b,  path1=path2, fileName1=paste(fileName2, "_boxPlot_NoFacet3",   sep="",  collapse=NULL),  height1=height2, width1=width2)

  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis ) ) +  
    geom_boxplot( outlier.size=0.1, notch=FALSE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot_NoFacet_withOutlier",   sep="",  collapse=NULL),  height1=height2, width1=width2-1)                             
 
  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleType  ) ) +  
    geom_boxplot( outlier.size=0.1, notch=FALSE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot_NoFacet_withOutlier2",   sep="",  collapse=NULL),  height1=height2, width1=width2)

  FigureTemp1 <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleType  ) ) +  
    geom_boxplot( outlier.size=0.1, notch=FALSE,  notchwidth = 0.1,  alpha=1  ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp1,  path1=path2, fileName1=paste(fileName2, "_boxPlot_NoFacet_withOutlier3",   sep="",  collapse=NULL),  height1=height2, width1=width2)
 
 
  FigureTemp4a <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis ) ) +  
    geom_violin(  colour = "red", fill="red"  ) + 
    geom_boxplot( outlier.shape=NA, outlier.size=0, size=0.6,    width=0.3,     alpha=0.00001, position=position_dodge(width=0.9)    ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4a,  path1=path2, fileName1=paste(fileName2, "_violinBoxPlot_NoFacet",   sep="",  collapse=NULL),  height1=height2, width1=width2-1)
   
  FigureTemp4a <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleType ) ) +  
    geom_violin(  colour = NA  ) + 
    geom_boxplot( outlier.shape=NA, outlier.size=0, size=0.6,    width=0.3,     alpha=0.00001, position=position_dodge(width=0.9)    ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4a,  path1=path2, fileName1=paste(fileName2, "_violinBoxPlot_NoFacet2",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  
  FigureTemp4a <- ggplot( DataFrame_Local, aes(x=sampleType, y=yAxis, fill=sampleType ) ) +  
    geom_violin(  colour = NA  ) + 
    geom_boxplot( outlier.shape=NA, outlier.size=0, size=0.6,    width=0.3,     alpha=0.00001, position=position_dodge(width=0.9)    ) +   
    stat_summary( position=position_dodge(width=0.9), fun.y=mean,  color="yellow4",  geom="point", shape=19, size=1.5, show.legend = FALSE) + 
    scale_fill_manual( values = colours2 ) +
    xlab(xLab2 ) + ylab( yLab2 ) + ggtitle( title2 )  + MyTheme_1_g(textSize1=12, hjust1=1, vjust1=1,  angle1=30 ) + ylim(Ymin2, Ymax2 )
  MySaveGgplot2_1_g(ggplot2Figure1=FigureTemp4a,  path1=path2, fileName1=paste(fileName2, "_violinBoxPlot_NoFacet3",   sep="",  collapse=NULL),  height1=height2, width1=width2)
  

}  



myOutPath_3_bed  = paste(outDir_g, "3-Formatted/BED",   sep="/")  

diff1_hyper_matrix1 = read.table(file= paste(myOutPath_3_bed, "diff1_hyper.bed",  sep="/" ), header = FALSE, sep = "\t", quote = "\"'", dec = ".")
diff1_hypo_matrix1  = read.table(file= paste(myOutPath_3_bed, "diff1_hypo.bed",   sep="/" ), header = FALSE, sep = "\t", quote = "\"'", dec = ".")

dim(diff1_hyper_matrix1)
dim(diff1_hypo_matrix1)
diff1_hyper_matrix1 = as.matrix( diff1_hyper_matrix1)   
diff1_hypo_matrix1  = as.matrix( diff1_hypo_matrix1)
dim(diff1_hyper_matrix1)
dim(diff1_hypo_matrix1)


diff1_hyper_length = as.numeric(diff1_hyper_matrix1[,3]) - as.numeric(diff1_hyper_matrix1[,2]) 
diff1_hypo_length  = as.numeric(diff1_hypo_matrix1[,3])  - as.numeric(diff1_hypo_matrix1[,2])  
length(diff1_hyper_length)
length(diff1_hypo_length)


vector_1 = c(diff1_hyper_length, diff1_hypo_length)
vector_2 = c( rep("hyper" , times=length(diff1_hyper_length)),   rep("hypo" , times=length(diff1_hypo_length)) )
matrix2 = cbind(vector_1, vector_2)
dim(matrix2) 


myOutPath_4  = paste(outDir_g, "4-LengthDis",   sep="/")  
if( ! file.exists(myOutPath_4)   ) { dir.create(myOutPath_4, recursive = TRUE) }

write.table(x=matrix2, file = paste(myOutPath_4, "/", "matrix", ".txt", sep=""), 
           append = FALSE, quote = FALSE, sep = "\t",
           eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE )

sink( file=paste(myOutPath_4, "/", "summary", ".txt", sep="") )
summary(vector_1)
sink() 

MyBoxViolinPlot_1_f(vector2=vector_1,   sampleType2=vector_2,  
                    colours2=c("red3", "cyan3"),   
                    path2=myOutPath_4,   fileName2="plots-1",  
                    title2="plots-1",  xLab2="DMRs",  yLab2="Length of DMRs",    
                    height2=4,   width2=3,   Ymin2=100, Ymax2=250)   

MyBoxViolinPlot_1_f(vector2=vector_1,   sampleType2=vector_2,  
                    colours2=c("red3", "cyan3"),   
                    path2=myOutPath_4,   fileName2="plots-2",  
                    title2="plots-2",  xLab2="DMRs",  yLab2="Length of DMRs",    
                    height2=4,   width2=3,   Ymin2=100, Ymax2=300)   


MyBoxViolinPlot_1_f(vector2=vector_1,   sampleType2=vector_2,  
                    colours2=c("red3", "cyan3"),   
                    path2=myOutPath_4,   fileName2="plots-3",  
                    title2="plots-3",  xLab2="DMRs",  yLab2="Length of DMRs",    
                    height2=4,   width2=3,   Ymin2=100, Ymax2=500)   



#############################################





## Part 7: Annotation of genomic regions. 
#############################################
myOutPath_3_chipseeker  = paste(outDir_g, "3-Formatted/forChIPseeker",   sep="/")  

myOutPath_5  = paste(outDir_g, "5-Annotation",   sep="/")  
if( ! file.exists(myOutPath_5)   ) { dir.create(myOutPath_5, recursive = TRUE) }  

myCommands_5 =  paste(  "Rscript  useChIPseeker.R  ", myOutPath_3_chipseeker,  "    " , myOutPath_5,  "     ", RefGenome_g ) 
system( myCommands_5 )
#############################################










