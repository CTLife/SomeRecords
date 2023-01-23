suppressPackageStartupMessages( library(ComplexHeatmap) ) 
suppressPackageStartupMessages( library(corrplot) )
suppressPackageStartupMessages( library(gplots) )
suppressPackageStartupMessages( library(stringr) )
suppressPackageStartupMessages( library(ggplot2) )


reset_outliers2 <- function(x, na.rm = TRUE ) {
  qnt <- quantile(x, probs=c(0.05, 0.95) , type=1,  na.rm = na.rm )  
  y <- x
  y[x < qnt[1] ] <- qnt[1]
  y[x > qnt[2] ] <- qnt[2]    
  y
}


myScaleMatrix2 <- function( matrix_temp8, upper_temp8 = 1, lower_temp8 = -1 ) {
    rawMatrix_2 = reset_outliers2(matrix_temp8)  
    rawMatrix_2 = lower_temp8 + (upper_temp8 - lower_temp8) * ( rawMatrix_2 - min(rawMatrix_2) )/( max(rawMatrix_2)- min(rawMatrix_2) )
    return(rawMatrix_2)
}


MyHeatmaps_1_cluster <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30, is.corr2=TRUE,  my_col2 ) {  ## perform clustering                                                           
  ##matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( corrplot(matrix2, method = "circle", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",   tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "square", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",   tl.col = "black", col = my_col2 )  )                                       
      print( corrplot(matrix2, method = "square", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "complete",  tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "number", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",   tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",   tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "complete",  tl.col = "black", col = my_col2 )  )   
  dev.off()
}  


MyHeatmaps_1_original <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30, is.corr2=TRUE,  my_col2 ) { ## original order                                                            
  ## matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( corrplot(matrix2, method = "circle", type = "full",   title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "square", type = "full",   title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", col = my_col2 )  )                                       
      print( corrplot(matrix2, method = "square", type = "lower",  title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "number", type = "full",   title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",   title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "lower",  title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", col = my_col2 )  )   
  dev.off()
}  


MyHeatmaps_2_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "red"),  bias = 2  )(20),       symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("green4", "white", "purple"),  bias = 2  )(20),  symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "red"),  bias = 2  )(20),       symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "purple"),  bias = 2  )(20),    symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "purple"),  bias = 2  )(20),    symm = FALSE,  scale =  "none"     )  )
  dev.off()
}  


MyHeatmaps_3_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )          
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20)   , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20)   , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20)        , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
      print( ComplexHeatmap::Heatmap(matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20)     , heatmap_legend_param = list(legend_height = unit(10, "cm")) )  )
  dev.off()
}  





#################################################
rawMatrix_1 <- read.table("Intervene_pairwise_frac_matrix.txt", header=TRUE, sep="\t", quote = "", comment.char = "") 
rownames( rawMatrix_1 ) = rawMatrix_1[,1]
rawMatrix_1 = rawMatrix_1[,-1]
rownames(rawMatrix_1)
colnames(rawMatrix_1)
dim(rawMatrix_1)

rawMatrix_1 = as.matrix(rawMatrix_1)


for(i in c(1:nrow(rawMatrix_1)) ) {
  rawMatrix_1[i,i] = 0
}
rawMatrix_1 = rawMatrix_1 * 100
rawMatrix_1 = round(rawMatrix_1, digits = 2)

outPath = "figures"
if( ! file.exists(outPath)  ) { dir.create(outPath,  recursive = TRUE)  }

my_col1=colorRampPalette( c(   "green4", "green", "yellow",  "yellow4", "purple",  "purple4"),  bias = 2  )
MyHeatmaps_1_original(matrix2=rawMatrix_1,  path2=outPath,   fileName2="1.heatmap.pdf",  height2=7,   width2=8, is.corr2=FALSE,   my_col2=my_col1(1000)  )

MyHeatmaps_2_f(matrix2=rawMatrix_1,  path2=outPath,   fileName2="2.heatmap.pdf",   height2=7,   width2=8 ) 


my_col2=colorRampPalette( c(    "grey", "yellow",  "yellow4", "purple",  "purple4") , bias = 2)  
pdf( file=paste(outPath, "3.heatmap.pdf", sep="/"),  height=8,   width=8 )
heatmap.2 (x= rawMatrix_1 ,
           # dendrogram control
           dendrogram = "none", 
           Rowv = FALSE,
           Colv="Rowv" , 
           symm = FALSE,
           # data scaling
           scale =  "none" ,
           na.rm=TRUE,
           # colors
           col=my_col2(40),
           trace = "none", 
           # cell labeling
           cellnote= as.matrix(rawMatrix_1) ,
           notecex=0.5,
           notecol="white",
           na.color=par("bg") 
)
dev.off() 


my_col3=colorRampPalette( c( "cyan4",  "cyan", "red",   "red4"),bias = 2 )
pdf( file = paste(outPath, "4.heatmap.pdf", sep="/"),  width=7, height=8  )
print( corrplot(rawMatrix_1, method = "square", type = "full",  title = "", is.corr = FALSE,  order = "original",  tl.col = "black", tl.srt = 90, col = my_col3(1000) )  )  
dev.off() 





my_col3=colorRampPalette( c( "grey",  "pink", "red",   "red4"),bias = 2 )
pdf( file = paste(outPath, "5.heatmap.pdf", sep="/"),  width=8, height=8  )
print( corrplot(rawMatrix_1, method = "square", type = "full",  title = "", is.corr = FALSE,  order = "original",  tl.col = "black", tl.srt = 90, col = my_col3(1000) )  )  
dev.off() 




