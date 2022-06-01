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
rawMatrix_1 <- read.table("numbers.txt", header=FALSE, sep="\t", quote = "", comment.char = "") 
dim(rawMatrix_1)
tail(rawMatrix_1)

rawMatrix_1A = rawMatrix_1 



outDir_g = "z.figures"
if( ! file.exists(outDir_g)   ) { dir.create(outDir_g, recursive = TRUE) }


p1 <- ggplot(rawMatrix_1A, aes(x=V2, y=log2(V1), fill=as.factor(V3) )) + geom_bar(stat="identity", color="black", position=position_dodge( )) +
        labs(x = "samples", y="Number mutation sites, log2") + scale_y_continuous(limits = c(0, 22), breaks = seq(0, 22, 2) )

MySaveGgplot2_1_g(ggplot2Figure1=p1, path1=outDir_g, fileName1="barPlot", height1=4, width1=7)



