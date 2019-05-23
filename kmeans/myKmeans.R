##################################################################################################################
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  myKmeans.R     2.Arithmetic_mean_matrix.txt   5     
         
args_g <- commandArgs(TRUE)
print("##########################")
print("args: ")
print(args_g[1])   
print(args_g[2])  
print("##########################")

input_g     = args_g[1];     ## Input matrix file
numClusters = args_g[2];     ## Number of clusters
numClusters = as.numeric(numClusters)
print( paste("numClusters:", numClusters, sep="") )  

# input_g =  "2.Arithmetic_mean_matrix.txt"
# numClusters = 5 

  
outDir_g = paste("Kmeans_k=", numClusters, sep="")
if( ! file.exists(input_g)    ) { print("##### Error-1 #####") }
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }
##################################################################################################################


 


###################
rawMatrix_1 <- read.table(input_g, header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_1)
colnames(rawMatrix_1)


rawMatrix_100 <- read.table("3B_mat-tiles.txt", header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_100)
colnames(rawMatrix_100)

rawMatrix_200 <- read.table("3C_regions-tiles.txt", header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_200)
colnames(rawMatrix_200)


nrow(rawMatrix_1) == nrow(rawMatrix_100)
nrow(rawMatrix_1) == nrow(rawMatrix_200)


#rownames(rawMatrix_1)
my_row_names = rawMatrix_1[,1]
  
rawMatrix_2 = rawMatrix_1[,-1]
rownames(rawMatrix_2) = my_row_names
colnames(rawMatrix_2)
#rownames(rawMatrix_2)
dim(rawMatrix_2)




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

rawMatrix_3 = rawMatrix_2 
for(i in c(1:nrow(rawMatrix_2)) ) {
  rawMatrix_3[i,] = myScaleMatrix2( matrix_temp8 = rawMatrix_2[i,], upper_temp8 = 1, lower_temp8 = -1  )
}

dim(rawMatrix_1)
dim(rawMatrix_2)
dim(rawMatrix_3)

write.table(x=rawMatrix_3, file = ( file = paste(outDir_g, "1.Row-scaled-matrix.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )



################
results_1 = kmeans(x=rawMatrix_3, centers=numClusters, iter.max = 1000, nstart =25)

sink(  paste(outDir_g, "2.kmeans_output.txt", sep="/") )
    print(results_1)
sink() 


results_1_cluster = results_1$cluster
write.table(x=results_1_cluster, file = ( file = paste(outDir_g, "3.cluster-of-each-row.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )

Num_each_cluster = rep(0, times=numClusters )
names_2 = rep(0, times=numClusters )
for(i in c(1:numClusters) ) {
    Num_each_cluster[i] = length( results_1_cluster[results_1_cluster==i] )
    names_2[i] = paste("k=", i, sep="")
}
names(Num_each_cluster) = names_2
write.table(x=Num_each_cluster, file = ( file = paste(outDir_g, "4.Number-rows-of-each-cluster.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )



results_1_cluster_index = order(results_1_cluster)

write.table(x=results_1_cluster_index, file = ( file = paste(outDir_g, "5.ordered_cluster_index.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


rawMatrix_4 = rawMatrix_3[results_1_cluster_index, ]

write.table(x=cbind(rawMatrix_4,results_1_cluster_index), file = ( file = paste(outDir_g, "6.Matrix_ordered.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


results_1_cluster4 = results_1_cluster[results_1_cluster_index]
write.table(x=results_1_cluster4, file = ( file = paste(outDir_g, "6.results_1_cluster4.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


my_col1=colorRampPalette( c( "blue",  "white", "red"  ), bias = 1 )

library(gplots)
pdf( file=paste(outDir_g, "7.heatmap.pdf",sep="/"),  height=5,   width=5 )
heatmap.2(x=as.matrix(rawMatrix_4),
           
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

for( i  in  c(1:numClusters) ) {
  my_temp1 = rawMatrix_200[ which(results_1_cluster==i), ]
  my_temp2 = my_row_names[ which(results_1_cluster==i)]
  my_temp3 = cbind( my_temp1[,-4],  my_temp2)  
  write.table(x=my_temp3, file = paste(for_regions, "/",  "cluster-", i, ".bed", sep="") , 
              append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = FALSE )
}



############
colnames(rawMatrix_100)
colnames(rawMatrix_4)

rawMatrix_100_X1NC         = rawMatrix_100[, 1:34 ]
rawMatrix_100_X2IVFfresh   = rawMatrix_100[, 35:60 ]
rawMatrix_100_X4ICSIfresh  = rawMatrix_100[, 61:90 ]
rawMatrix_100_X3IVFfrozen  = rawMatrix_100[, 91:118 ]
rawMatrix_100_X5ICSIfrozen = rawMatrix_100[, 119:137 ]

my_matrix1 = matrix(data = NA, nrow = nrow(rawMatrix_100), ncol = 1)
for( i  in  c(1:ncol(rawMatrix_4)) ) {
     bool_1 = grepl(pattern=colnames(rawMatrix_4)[i], x="rawMatrix_100_X1NC")
     bool_2 = grepl(pattern=colnames(rawMatrix_4)[i], x="rawMatrix_100_X2IVFfresh")
     bool_3 = grepl(pattern=colnames(rawMatrix_4)[i], x="rawMatrix_100_X4ICSIfresh")
     bool_4 = grepl(pattern=colnames(rawMatrix_4)[i], x="rawMatrix_100_X3IVFfrozen")
     bool_5 = grepl(pattern=colnames(rawMatrix_4)[i], x="rawMatrix_100_X5ICSIfrozen")
     if(bool_1) { my_matrix1 = cbind(my_matrix1, rawMatrix_100_X1NC)          }
     if(bool_2) { my_matrix1 = cbind(my_matrix1, rawMatrix_100_X2IVFfresh)    }
     if(bool_4) { my_matrix1 = cbind(my_matrix1, rawMatrix_100_X3IVFfrozen)   }
     if(bool_3) { my_matrix1 = cbind(my_matrix1, rawMatrix_100_X4ICSIfresh)   }
     if(bool_5) { my_matrix1 = cbind(my_matrix1, rawMatrix_100_X5ICSIfrozen)  }
}
my_matrix1 = my_matrix1[,-1]
dim(my_matrix1)

for_MeLevel = paste(outDir_g, "MeLevel_eachNeonate", sep="/")
if( ! file.exists(for_MeLevel)    ) { dir.create(for_MeLevel, recursive = TRUE) }

my_matrix2 = my_matrix1[results_1_cluster_index, ]
write.table(x=my_matrix2, file = paste(for_MeLevel, "matrix_sorted.txt",   sep="/") , 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE )

my_matrix3 = my_matrix2
for(i in c(1:nrow(my_matrix2)) ) {
  my_matrix3[i,] = myScaleMatrix2( matrix_temp8 = my_matrix2[i,], upper_temp8 = 10, lower_temp8 = -10  )
}
write.table(x=my_matrix3, file = paste(for_MeLevel, "matrix_sorted_scaled.txt",   sep="/") , 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE )


my_col1=colorRampPalette( c( "blue",  "white", "white", "white", "red"  ), bias = 1 )

pdf( file=paste(for_MeLevel, "heatmap.pdf",sep="/"),  height=5,   width=15 )
heatmap.2(x=as.matrix(my_matrix3),
          
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











##############################################
results_1_cluster_index2 = c(    which(results_1_cluster==1), which(results_1_cluster==5), 
                                 which(results_1_cluster==6), which(results_1_cluster==2), 
                                 which(results_1_cluster==3), which(results_1_cluster==4),
                                 which(results_1_cluster==7), which(results_1_cluster==9),  
                                 which(results_1_cluster==10), which(results_1_cluster==8), 
                                 which(results_1_cluster==12), which(results_1_cluster==11) )

length(results_1_cluster_index)
length(results_1_cluster_index2)


write.table(x=results_1_cluster_index2, file = ( file = paste(outDir_g, "105.ordered_cluster_index2.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


rawMatrix_42 =rawMatrix_3[results_1_cluster_index2, ]
write.table(x=cbind(rawMatrix_42,results_1_cluster_index2), file = ( file = paste(outDir_g, "106.Matrix_ordered2.txt",sep="/") ), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = TRUE,  col.names = TRUE )


my_col1=colorRampPalette( c( "blue",  "white", "red"  ), bias = 1 )

library(gplots)
pdf( file=paste(outDir_g, "107.heatmap.pdf",sep="/"),  height=5,   width=5 )
heatmap.2(x=as.matrix(rawMatrix_42),
          
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




