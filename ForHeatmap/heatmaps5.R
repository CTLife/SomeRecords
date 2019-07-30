args <- commandArgs(TRUE)

print("args: ")
print(args[1])   
print("#############")
file1 = args[1]; 


## file1 = "1.GO-BP.DMRs.IVF-ET_vs_CTRL.txt"

suppressPackageStartupMessages( library(ComplexHeatmap) ) 
suppressPackageStartupMessages( library(corrplot) )
library(gplots)

reset_outliers2 <- function(x, na.rm = TRUE ) {
  qnt <- quantile(x, probs=c(0.1, 0.9) , type=1,  na.rm = na.rm )  
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



MyHeatmaps_1_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30, is.corr2=TRUE,  my_col2 ) {                                                             
  matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( corrplot(matrix2, method = "circle", type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "circle", type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )                                       
      print( corrplot(matrix2, method = "number", type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "pie",    type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "upper", title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "hclust",  hclust.method = "ward.D2",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  dev.off()
  
}  



MyHeatmaps_1A_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30, is.corr2=TRUE,  my_col2 ) {                                                             
  ## matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( corrplot(matrix2, method = "circle", type = "full",  title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "circle", type = "upper", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )                                       
      print( corrplot(matrix2, method = "number", type = "full", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "pie",    type = "full", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full", title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
      print( corrplot(matrix2, method = "color",  type = "full",  title = "", is.corr = is.corr2,  order = "original",  tl.col = "black", tl.srt = 45, col = my_col2 )  )   
  dev.off()
  
}  



MyHeatmaps_2_f <- function(matrix2,  path2,   fileName2,   height2=30,   width2=30 ) { 
  matrix2 = myScaleMatrix2( matrix2 )  
  pdf( file = paste(path2, fileName2, sep="/"),  width=width2, height=height2  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "red"))(20),       symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("green4", "white", "purple"))(20),  symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "red"))(20),       symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("cyan", "white", "purple"))(20),    symm = FALSE,  scale =  "none"     )  )
      print( heatmap(x = matrix2,  col = colorRampPalette(c("blue", "white", "purple"))(20),    symm = FALSE,  scale =  "none"     )  )
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
rawMatrix_1 <- read.table(file1, header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_1)

rawMatrix_1[ is.na(rawMatrix_1) ] = 1
rawMatrix_2 = -log10(rawMatrix_1)


outPath = paste("Figures", file1, sep=".")
if( ! file.exists(outPath) ) { dir.create(outPath, recursive = TRUE) }

my_col1=colorRampPalette( c("grey", "green4",  "yellow",  "yellow" ), bias = 1.5 )

pdf( file= paste(outPath, "heatmap-3.pdf", sep="/"),  height=5,   width=5 )
heatmap.2 (x=as.matrix(rawMatrix_2),
           
           # dendrogram control
           dendrogram = "none", 
           Rowv = FALSE,
           Colv="Rowv" , 
           symm = FALSE,
           
           # data scaling
           scale =  "none" ,
           na.rm=TRUE,
           
           # colors
           col=my_col1(20),
           trace = "none", 
           
           # cell labeling
           #cellnote= as.matrix(rawMatrix_2) ,
           notecex=1.0,
           notecol="white",
           na.color=par("bg") 
            
)
dev.off() 







