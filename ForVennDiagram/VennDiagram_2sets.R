library(VennDiagram)


spe1 <-  732 
common <- 25093
spe2 <- 4335


x1 <- seq(from = 1, to = (spe1+common), by = 1)
length(x1)

x2 <- seq(from = spe1+1, to = (spe1+common+spe2), by = 1)
length(x2)


myFile <- "group1_vs_group2.svg"

venn.diagram(x=list("group1"=x1,   "group2"=x2), 
             filename=myFile, 
             height = 2,   width = 2,  
             resolution =500,   imagetype = "svg",   units = "px", 
             compression ="lzw",  na = "stop",  
             main = "group1 vs group2", sub = NULL, 
             main.pos= c(0.5, 1.05),  main.fontface = "plain",
             main.fontfamily = "serif",  main.col ="black",
             main.cex = 1,  main.just = c(0.5, 1), 
             sub.pos = c(0.5, 1.05),  sub.fontface = "plain", 
             sub.fontfamily ="serif", sub.col = "black", sub.cex = 1,
             sub.just =c(0.5, 1), #category.names = names(x), 
             force.unique =TRUE, print.mode = "raw", sigdigs = 3, 
             direct.area =FALSE, area.vector = 0, hyper.test = TRUE,
             fill=c("green3", "yellow3"), cat.col=c("green3", "yellow3"), lwd=0, margin=0.1 )





