matrix_1 <- read.table("1-raw.nsmb.2660-S2",  header=FALSE,  sep="\t",  quote = "",  comment.char = "") 
matrix_1[1:10,1:10]
dim(matrix_1)

matrix_2 <- matrix_1[-1,]
matrix_2[1:10,1:10]
dim(matrix_2)

num_row    <- nrow(matrix_2)
num_column <- ncol(matrix_2)

rownames(matrix_1) <- c()
colnames(matrix_1) <- c()

column_names <- as.vector( unlist(matrix_1[1,]) )
length(column_names)
dim(column_names)

column_names <- as.vector( column_names[-1] )
length(column_names)
dim(column_names)



matrix_3 <- matrix(data = NA, nrow = 0, ncol = 3, byrow = FALSE, dimnames = NULL)
dim(matrix_3)

rep(x=column_names[1], times = num_row, length.out = NA, each = 1) 



for (i in 2:num_column) {  ## include start and end.
  gene1    <- as.character( matrix_2[,1] )
  sample2  <- rep(x=column_names[i-1], times = num_row, length.out = NA, each = 1) 
  value3   <- as.character( matrix_2[,i] )
  length(gene1)
  length(sample2)
  length(value3)
  matrix_temp <- cbind(gene1, sample2, value3)
  matrix_3 <- rbind(matrix_3, matrix_temp)
}




dim(matrix_3)
num_row*(num_column-1)

matrix_3[1:100, ]
matrix_3[2506536, ]

write.table(x=matrix_3, file = "matrix_3.txt", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")







