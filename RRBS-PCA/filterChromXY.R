## example:  Rscript   filterChromXY.R     1-Coverage-CpG  1_A42C-girl-IVF-frozen_Rep1.bismark.cov   



args <- commandArgs(TRUE)
print("args: ")
print(args[1])   
print(args[2])   
print("#############")

inputDir   = args[1];
inputFile  = args[2];     ## input file
#  inputDir   = "1-Coverage-CpG";
#  inputFile  = "1_A42C-girl-IVF-frozen_Rep1.bismark.cov";     ## input file


outPath1   = "10A-rmXY"
outPath2   = "10B-onlyXY"
outPath3   = "10C-onlyX"
outPath4   = "10D-onlyY"
outPath5   = "10E-otherLog"
if( ! file.exists(outPath1) ) { dir.create(outPath1, recursive = TRUE)  }
if( ! file.exists(outPath2) ) { dir.create(outPath2, recursive = TRUE)  }       
if( ! file.exists(outPath3) ) { dir.create(outPath3, recursive = TRUE)  }
if( ! file.exists(outPath4) ) { dir.create(outPath4, recursive = TRUE)  }
if( ! file.exists(outPath5) ) { dir.create(outPath5, recursive = TRUE)  }


###########################################################################################
matrix_1 <- read.table(paste(inputDir, inputFile, sep="/") , header=FALSE, sep="\t",  quote = "",  comment.char = "") 
myChrom = matrix_1[,1]
matrix_1A <- matrix_1[(myChrom!="chrX" & myChrom!="chrY"), ] 
matrix_1B <- matrix_1[(myChrom=="chrX" | myChrom=="chrY"), ] 
matrix_1C <- matrix_1[myChrom=="chrX", ] 
matrix_1D <- matrix_1[myChrom=="chrY", ] 


write.table(matrix_1A, file = paste(outPath1, inputFile, sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)
write.table(matrix_1B, file = paste(outPath2, inputFile, sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)
write.table(matrix_1C, file = paste(outPath3, inputFile, sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)
write.table(matrix_1D, file = paste(outPath4, inputFile, sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names =  FALSE,   col.names =  FALSE)

sink( file = paste(outPath5, "/",  "1-dimensions-", inputFile, sep="")  )
print("raw:")
print( dim(matrix_1) )
print("rmXY:")
print( dim(matrix_1A) )
print("onlyXY:")
print( dim(matrix_1B) )
print("onlyX:")
print( dim(matrix_1C) )
print("onlyY:")
print( dim(matrix_1D) )
sink() 





