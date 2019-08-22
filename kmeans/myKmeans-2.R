##################################################################################################################
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  myKmeans.R     3B_mat-tiles.txt     3C_regions-tiles.txt     6     
    

library(factoextra) # clustering algorithms & visualization
library(flexclust)

args_g <- commandArgs(TRUE)
print("##########################")
print("args: ")
print(args_g[1])   
print(args_g[2])  
print(args_g[3])  
print("##########################")

input_matrix_g= args_g[1];     ## Input matrix file
input_region_g= args_g[2];     ## Input bed file
numClusters_g = args_g[3];     ## Number of clusters
numClusters_g = as.numeric(numClusters_g)
print( paste("numClusters:", numClusters_g, sep=" ") )  

# input_matrix_g = "3B_mat-tiles.txt"
# input_region_g = "3C_regions-tiles.txt"
# numClusters_g  = 3 

outDir_g = paste("Kmeans_k=", numClusters_g, sep="")
if( ! file.exists(input_matrix_g)    ) { print("##### Error-1 #####") }
if( ! file.exists(input_region_g)    ) { print("##### Error-2 #####") }
if( ! file.exists(outDir_g)          ) { dir.create(outDir_g, recursive = TRUE) }
##################################################################################################################


 


###################
rawMatrix_100 <- read.table(input_matrix_g, header=TRUE, sep="\t", quote = "", comment.char = "") 
rawMatrix_200 <- read.table(input_region_g, header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_100)
dim(rawMatrix_200)
head(rawMatrix_100)
head(rawMatrix_200)
nrow(rawMatrix_100) == nrow(rawMatrix_200)  
colNames_samples = colnames(rawMatrix_100)

NC_cols        = c(1:34)
IVFfresh_cols  = c(35:60)
ICSIfresh_cols = c(61:90)
IVFfrozen_cols = c(91:118)
ICSIfrozen_cols= c(119:137)
NC_cols 
IVFfresh_cols 
ICSIfresh_cols 
IVFfrozen_cols 
ICSIfrozen_cols 
length( NC_cols )
length( IVFfresh_cols )
length( ICSIfresh_cols )
length( IVFfrozen_cols )
length( ICSIfrozen_cols )

selected_cols = c(NC_cols, IVFfresh_cols, ICSIfresh_cols)
rawMatrix_101 = rawMatrix_100[, selected_cols]  
dim(rawMatrix_101)


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



## scale each row
rawMatrix_102 = rawMatrix_101
for(i in c(1:nrow(rawMatrix_101)) ) {
  ave_temp = mean( as.numeric(rawMatrix_101[i,]), na.rm=TRUE)
  temp_1 = rawMatrix_101[i,]
  temp_1[is.na(temp_1)] = ave_temp
  rawMatrix_102[i,] = myScaleMatrix2( matrix_temp8 = temp_1, upper_temp8 = 1, lower_temp8 = -1  )
}
dim(rawMatrix_100)
dim(rawMatrix_101)
dim(rawMatrix_102)


rawMatrix_103 = cbind(rawMatrix_200, rawMatrix_102)
write.table(x=rawMatrix_103, file = ( file = paste(outDir_g, "1.Row-scaled-matrix.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE )






################
## distCor, distJaccard
cluster_results_1 = flexclust::kcca(x=rawMatrix_102, k=3, family=kccaFamily(which=NULL, dist=distCor), 
                                    weights=NULL, group=NULL, control=NULL, simple=FALSE, save.data=FALSE)
sink(  paste(outDir_g, "2.kmeans_output.txt", sep="/") )
  print(cluster_results_1)
sink() 


results_1_cluster = cluster_results_1@cluster
write.table(x=results_1_cluster, file = ( file = paste(outDir_g, "3.cluster-of-each-row.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )

Num_each_cluster = rep(0, times=numClusters_g )
names_2 = rep(0, times=numClusters_g )
for(i in c(1:numClusters_g) ) {
  Num_each_cluster[i] = length( results_1_cluster[results_1_cluster==i] )
  names_2[i] = paste("k=", i, sep="")
}
names(Num_each_cluster) = names_2
write.table(x=Num_each_cluster, file = ( file = paste(outDir_g, "4.Number-rows-of-each-cluster.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )



results_1_cluster_index = order(results_1_cluster)
write.table(x=results_1_cluster_index, file = ( file = paste(outDir_g, "5.ordered_cluster_index.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


rawMatrix_104 = rawMatrix_102[results_1_cluster_index, ]
write.table(x=cbind(rawMatrix_104,results_1_cluster_index), file = ( file = paste(outDir_g, "6.Matrix_ordered.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


results_1_cluster4 = results_1_cluster[results_1_cluster_index]
write.table(x=results_1_cluster4, file = ( file = paste(outDir_g, "7.sorted.cluster-of-each-row.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


my_col1=colorRampPalette( c( "blue",  "white", "red"  ), bias = 1 )

library(gplots)

pdf( file=paste(outDir_g, "8.heatmap.pdf",sep="/"),  height=10,   width=5 )
heatmap.2(x=as.matrix(rawMatrix_104),
           
           # dendrogram control
           dendrogram = "none", 
           Rowv = FALSE,
           Colv="Rowv" , 
           symm = FALSE,
           
           # data scaling
           scale =  "none" ,
           na.rm=TRUE,
           
           # colors
           col=my_col1(40),
           trace = "none", 
           
           # cell labeling
           #cellnote= as.matrix(rawMatrix_4) ,
           notecex=1.0,
           notecol="white",
           na.color=par("bg") 
           
)
dev.off() 

png( file=paste(outDir_g, "8.heatmap.png",sep="/"),  height=1000,   width=500 )
heatmap.2(x=as.matrix(rawMatrix_104),
          
          # dendrogram control
          dendrogram = "none", 
          Rowv = FALSE,
          Colv="Rowv" , 
          symm = FALSE,
          
          # data scaling
          scale =  "none" ,
          na.rm=TRUE,
          
          # colors
          col=my_col1(40),
          trace = "none", 
          
          # cell labeling
          #cellnote= as.matrix(rawMatrix_4) ,
          notecex=1.0,
          notecol="white",
          na.color=par("bg") 
          
)
dev.off() 


for_regions = paste(outDir_g, "genomicRegions/1-bed", sep="/")
if( ! file.exists(for_regions)    ) { dir.create(for_regions, recursive = TRUE) }

for( i  in  c(1:numClusters_g) ) {
  my_temp1 = rawMatrix_200[ which(results_1_cluster==i), ]
  write.table(x=my_temp1, file = paste(for_regions, "/",  "cluster-", i, ".bed", sep="") , 
              append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE )
}









