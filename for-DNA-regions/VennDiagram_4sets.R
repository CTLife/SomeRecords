library(VennDiagram)


myA = 4844  # ICSI_diff1_hyper
myB = 5820   # ICSI_diff1_hypo
myC = 6317   # IVF_diff1_hyper
myD = 5696   # IVF_diff1_hypo

myBD = 851
myBC = 992
myAD = 810
myAC = 870


x_A <- paste('A',  1:myA, sep='')
x_B <- paste('B',  1:myB, sep='')
x_C <- paste('C',  1:myC, sep='')
x_D <- paste('D',  1:myD, sep='')
x_BD  <- paste('BD',  1:myBD, sep='')
x_BC  <- paste('BC',  1:myBC, sep='')
x_AD  <- paste('AD',  1:myAD, sep='')
x_AC  <- paste('AC',  1:myAC, sep='')

x1 <- c(x_A, x_AD, x_AC )
x2 <- c(x_B, x_BD, x_BC)
x3 <- c(x_C, x_BC, x_AC )
x4 <- c(x_D, x_BD, x_AD ) 

myFile <- "A_vs_B_vs_C_vs_D.svg"

venn.diagram(x=list("ICSI_hyper"=x1,   "ICSI_hypo"=x2,  "IVF_hyper"=x3,  "IVF_hypo"=x4), 
             filename=myFile, 
             height = 3,   width = 3,  
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
fill=c("skyblue", "blue3", "pink", "red"), cat.col=c("skyblue", "blue3", "pink", "red"), lwd=0, margin=0.1 )


