library(VennDiagram)


myA = 18245
myB = 403
myC = 100
myAB = 1097
myAC = 70
myBC = 28
myABC = 263


x_A <- paste('A',  1:myA, sep='')
x_B <- paste('B',  1:myB, sep='')
x_C <- paste('C',  1:myC, sep='')
x_AB  <- paste('AB',  1:myAB, sep='')
x_AC  <- paste('AC',  1:myAC, sep='')
x_BC  <- paste('BC',  1:myBC, sep='')
x_ABC <- paste('ABC',  1:myABC, sep='')

x1 <- c(x_A, x_AB, x_AC, x_ABC)
x2 <- c(x_B, x_AB, x_BC, x_ABC)
x3 <- c(x_C, x_BC, x_AC, x_ABC)


myFile <- "A_vs_B_vs_C.svg"

venn.diagram(x=list("group1"=x1,   "group2"=x2,  "group3"=x3), 
             filename=myFile, 
             height = 2,   width = 2,  
             resolution =500,   imagetype = "svg",   units = "px", 
             compression ="lzw",  na = "stop",  
             main = "A_vs_B_vs_C", sub = NULL, 
             main.pos= c(0.5, 1.05),  main.fontface = "plain",
             main.fontfamily = "serif",  main.col ="black",
             main.cex = 1,  main.just = c(0.5, 1), 
             sub.pos = c(0.5, 1.05),  sub.fontface = "plain", 
             sub.fontfamily ="serif", sub.col = "black", sub.cex = 1,
             sub.just =c(0.5, 1), #category.names = names(x), 
             force.unique =TRUE, print.mode = "raw", sigdigs = 3, 
             direct.area =FALSE, area.vector = 0, hyper.test = TRUE,
fill=c("green3", "yellow3", "blue"), cat.col=c("green3", "yellow3", "blue"), lwd=0, margin=0.1 )





