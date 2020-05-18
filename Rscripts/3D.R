library("car")
library("plot3D")
setwd("/media/yp/YP-FOREVER1/SAC-seq_Downstream/Figures/1")


matrix1 = read.table(file="formatted.txt", header = TRUE, sep = "\t", quote = "\"'" )
myregions =  matrix1$type2
myBool = (myregions<3)


matrix1A = matrix1[myBool, ]
matrix1B = matrix1[!myBool, ]


##my_col1=colorRampPalette( c( "green4",    "yellow" ), bias = 1.5 )


pdf("1.pdf" )
with(matrix1, scatter3D(x = gene, y = type2, z = day, colvar = mutation,
	pch = 16, cex = 1.5, xlab = "gene", ylab = "region", 	zlab = "day",  
	d = 2, type = "h",  theta=70, phi=0, bty = "g",  
	col = ramp.col(c("black", "yellow", "yellow", "cyan",  "cyan",   "pink", "pink", "red", "red", "red","red", "red4", "red4","red4", "red4")), 
	colkey = list(length = 0.5, width = 0.5, cex.clab = 1 )  )
)
dev.off()





pdf("2.pdf" )
with(matrix1A, scatter3D(x = gene, y = type2, z = day, colvar = mutation,
	pch = 16, cex = 1.5, xlab = "gene", ylab = "region", 	zlab = "day",  
	d = 2, type = "h",  theta=50, phi=0, bty = "g",  
	col = ramp.col(c("black", "yellow", "yellow", "cyan",  "cyan",   "pink", "pink", "red", "red", "red","red", "red4", "red4","red4", "red4")), 
	colkey = list(length = 0.5, width = 0.5, cex.clab = 1 )  )
)
with(matrix1B, scatter3D(x = gene, y = type2, z = day, colvar = mutation,
	pch = 16, cex = 1.5, xlab = "gene", ylab = "region", 	zlab = "day",  
	d = 2, type = "h",  theta=50, phi=0, bty = "g",  
	col = ramp.col(c("black", "yellow", "yellow", "cyan",  "cyan",   "pink", "pink", "red", "red", "red","red", "red4", "red4","red4", "red4")), 
	colkey = list(length = 0.5, width = 0.5, cex.clab = 1 )  )
)
dev.off()


