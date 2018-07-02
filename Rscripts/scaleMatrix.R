##################################################################################################################
## Suffixes of all self-defined global variables must be "_g".
## Example:  
## Rscript  scaleMatrix.R     10A_pearsonCorrelation-regions.txt     
         
args_g <- commandArgs(TRUE)
print("##########################")
print("args: ")
print(args_g[1])   
print("##########################")
input_g = args_g[1];     ## input matrix file
print(input_g) 

  
outDir_g = paste(input_g, ".scaled.txt", sep="") 
##################################################################################################################


 
 





###################
rawMatrix_1 <- read.table(input_g, header=TRUE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_1)


rawMatrix_2 = rawMatrix_1[,-c(1,2)]
dim(rawMatrix_2)


upper =  1
lower = -1
rawMatrix_2 = lower + (upper - lower) * ( rawMatrix_2 - min(rawMatrix_2) )/( max(rawMatrix_2)- min(rawMatrix_2) )
rawMatrix_3 = cbind(rawMatrix_1[,c(1,2)],  rawMatrix_2)

write.table(x=rawMatrix_3, file = outDir_g,  append = FALSE, quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = TRUE, qmethod = c("escape", "double"),
                 fileEncoding = "")








