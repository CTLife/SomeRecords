suppressPackageStartupMessages( library(optparse)  )  ## To run the script in command lines.
suppressPackageStartupMessages( library(tidyverse) )  ## ggplot2 and others.
suppressPackageStartupMessages( library(ggplot2) )  ## ggplot2 and others.

getParameters_f <- function() {
  option_list_Local <- list(   ## Options list with associated default value.  16 options.
    optparse::make_option(opt_str=c("-i", "--input"),  
                          default="GY1", 
                          type="character",   dest="input",
                          help="Name of the input file. [default: %default]." )                                            
  )
  
  ## Now parse the command line to check which option is given and get associated values.
  parser_Local <- optparse::OptionParser(usage="usage: %prog [options]",
                                         option_list=option_list_Local, 
                                         description="For the output of VarScan.",                             
                                         epilogue="For comments, bug reports etc..., please visit https://github.com/CTLife/PipiDABS or contact Yong Peng <yongp@outlook.com>."
  )
  opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
  return(opt_Local)
}
##############################################################################################################################################################################################





##############################################################################################################################################################################################
opt_g = getParameters_f()  

inputFile_g  <- opt_g$input
print("########################")
print(inputFile_g)
print("########################")

rm(getParameters_f)  
rm(opt_g)  

options(digits=10)
continue_on_error_g <- function() {
  print( "NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'. " )
}
options( error=continue_on_error_g )  ## This option is very important.

## inputFile_g = "GY1"
outDir_g = paste("Results", inputFile_g, sep="_" )
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g,   recursive = TRUE) }
##############################################################################################################################################################################################

MySaveGgplot2_1_g <- function(ggplot2Figure1,  path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggplot2::ggsave(filename=paste(SVG1, "/", fileName1, ".svg", sep=""),  plot = last_plot(), device = "svg",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)
  ggplot2::ggsave(filename=paste(PDF1, "/", fileName1, ".pdf", sep=""),  plot = last_plot(), device = "pdf",   path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)
  ggplot2::ggsave(filename=paste(EPS1, "/", fileName1, ".eps", sep=""),  plot = last_plot(), device =cairo_ps, path = NULL, scale = 1, width = width1, height = height1, units = "in", dpi = 3000, limitsize = FALSE)
}







DF_1 <- read.table( paste(inputFile_g, "A2C.bed", sep="/"), header=F,   sep="\t" )  
DF_2 <- read.table( paste(inputFile_g, "A2G.bed", sep="/"), header=F,   sep="\t" )  
DF_3 <- read.table( paste(inputFile_g, "A2T.bed", sep="/"), header=F,   sep="\t" )  
DF_4 <- read.table( paste(inputFile_g, "C2A.bed", sep="/"), header=F,   sep="\t" )  
DF_5 <- read.table( paste(inputFile_g, "C2G.bed", sep="/"), header=F,   sep="\t" )  
DF_6 <- read.table( paste(inputFile_g, "C2T.bed", sep="/"), header=F,   sep="\t" )  

dim(DF_1)
dim(DF_2)
dim(DF_3)
dim(DF_4)
dim(DF_5)
dim(DF_6)

distance_1 = as.numeric( DF_1[,21] )
distance_2 = as.numeric( DF_2[,21] )
distance_3 = as.numeric( DF_3[,21] )
distance_4 = as.numeric( DF_4[,21] )
distance_5 = as.numeric( DF_5[,21] )
distance_6 = as.numeric( DF_6[,21] )

distance = c(distance_1, distance_2, distance_3, distance_4, distance_5, distance_6)
types = c(  rep("A_to_C", length(distance_1)) , rep("A_to_G", length(distance_2)) ,rep("A_to_T", length(distance_3)) ,
            rep("C_to_A", length(distance_4)) , rep("C_to_G", length(distance_5)) ,rep("C_to_T", length(distance_6))  )

distance[distance>2000] = 2000
dataframeB = data.frame( xAxis=distance,  sampleType=types )

p1 <- ggplot(data=dataframeB, aes(x=xAxis,  colour=sampleType) )   +  xlab("distance") + ylab("Probability density") +  
      geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE   ) +
      scale_colour_manual( values= c("red", "orange", "blue", "grey80", "grey40", "grey20")   ) +
      scale_x_continuous(breaks=  seq(from = 0, to = 2000, by =100 ) )
MySaveGgplot2_1_g(ggplot2Figure1=p1, path1=outDir_g , fileName1="1.density", height1=3.5, width1=7)




distance[distance>1000] =1000
dataframeB = data.frame( xAxis=distance,  sampleType=types )

p1 <- ggplot(data=dataframeB, aes(x=xAxis,  colour=sampleType) )   +  xlab("distance") + ylab("Probability density") +  
  geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE   ) +
  scale_colour_manual( values= c("red", "orange", "blue", "grey80", "grey40", "grey20")   ) +
  scale_x_continuous(breaks= seq(from = 0, to = 1000, by =100 ) )
MySaveGgplot2_1_g(ggplot2Figure1=p1, path1=outDir_g , fileName1="2.density", height1=3.5, width1=7)




distance[distance>800] = 800
dataframeB = data.frame( xAxis=distance,  sampleType=types )

p1 <- ggplot(data=dataframeB, aes(x=xAxis,  colour=sampleType) )   +  xlab("distance") + ylab("Probability density") +  
  geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE   ) +
  scale_colour_manual( values= c("red", "orange", "blue", "grey80", "grey40", "grey20")   ) +
  scale_x_continuous(breaks= seq(from = 0, to = 800, by =100 ) )
MySaveGgplot2_1_g(ggplot2Figure1=p1, path1=outDir_g , fileName1="3.density", height1=3.5, width1=7)




distance[distance>500] = 500
dataframeB = data.frame( xAxis=distance,  sampleType=types )

p1 <- ggplot(data=dataframeB, aes(x=xAxis,  colour=sampleType) )   +  xlab("distance") + ylab("Probability density") +  
  geom_density(mapping = NULL, data = NULL, stat = "density", position = "identity", na.rm = FALSE   ) +
  scale_colour_manual( values= c("red", "orange", "blue", "grey80", "grey40", "grey20")   ) +
  scale_x_continuous(breaks= seq(from = 0, to = 500, by =100 ) )
MySaveGgplot2_1_g(ggplot2Figure1=p1, path1=outDir_g , fileName1="4.density", height1=3.5, width1=7)





