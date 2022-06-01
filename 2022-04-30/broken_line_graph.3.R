library(ggplot2)

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


###################
rawMatrix_1 <- read.table("matrix.txt", header=FALSE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_1)
tail(rawMatrix_1)

rawMatrix_1A = rawMatrix_1 

bool2 = ( rawMatrix_1A[,2]==2 | rawMatrix_1A[,2]==4 | rawMatrix_1A[,2]==6  )
rawMatrix_2A = rawMatrix_1A[bool2,]
dim(rawMatrix_2A)

bool3 = ( rawMatrix_1A[,2]==1 | rawMatrix_1A[,2]==3 | rawMatrix_1A[,2]==5  )
rawMatrix_3A = rawMatrix_1A[bool3,]
dim(rawMatrix_3A)



outDir_g = "z.figures"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }


p1 <- ggplot(rawMatrix_1A, aes(x=V3, y=V1, group=as.factor(V2) )) + geom_line(aes(color=as.factor(V2)) )+ labs(x = "Number of factors", y="Number of Peaks") +
      geom_point(aes(color=as.factor(V2)))  + scale_color_manual(values=c("cyan", "red",  "green4", "yellow4",  "purple" , "blue")) +
      scale_x_continuous(breaks=c(0,  2, 4, 6, 8, 10, 12, 14, 16, 18, 20) )
MySaveGgplot2_1_g(ggplot2Figure1=p1, path1=outDir_g, fileName1="1.all6", height1=3.5, width1=7)


p2 <- ggplot(rawMatrix_2A, aes(x=V3, y=V1, group=as.factor(V2) )) + geom_line(aes(color=as.factor(V2)) )+ labs(x = "Number of factors", y="Number of Peaks") +
      geom_point(aes(color=as.factor(V2)))  + scale_color_manual(values=c("red", "yellow4", "blue"  )) +
      scale_x_continuous(breaks=c(0,  2, 4, 6, 8, 10, 12, 14, 16, 18, 20) )
MySaveGgplot2_1_g(ggplot2Figure1=p2, path1=outDir_g, fileName1="2.only3", height1=3.5, width1=7)


p3 <- ggplot(rawMatrix_3A, aes(x=V3, y=V1, group=as.factor(V2) )) + geom_line(aes(color=as.factor(V2)) )+ labs(x = "Number of factors", y="Number of Peaks") +
      geom_point(aes(color=as.factor(V2)))  + scale_color_manual(values=c("red", "yellow4", "blue"  )) +
      scale_x_continuous(breaks=c(0,  2, 4, 6, 8, 10, 12, 14, 16, 18, 20) )
MySaveGgplot2_1_g(ggplot2Figure1=p3, path1=outDir_g, fileName1="3.only3", height1=3.5, width1=7)



 
