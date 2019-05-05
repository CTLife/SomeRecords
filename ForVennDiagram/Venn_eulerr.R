 library(eulerr) 
 library(Cairo) 


VennDiag1 <- euler( c("A" = 2875,  "C" = 1276,  "B" = 1202,   
                     "A&B" = 701, "B&C" = 190,  "A&C" = 342, "A&B&C" = 269), shape = "ellipse")

VennDiag2 <- euler( c("A" = 1475,  "C" = 1347,  "B" = 1576,   
                     "A&B" = 342, "B&C" = 178,  "A&C" = 154, "A&B&C" = 87), shape = "ellipse" )


pdf("1.pdf")
    plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("red3", "yellow3", "cyan3"), col=c("red3", "yellow3", "cyan3") )
    plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("purple3", "orange3", "green3"), col=c("purple3", "orange3", "green3") )
dev.off()



svg("1a.svg")
    plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("red3", "yellow3", "cyan3"), col=c("red3", "yellow3", "cyan3") )
dev.off()
svg("1b.svg")
    plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("purple3", "orange3", "green3"), col=c("purple3", "orange3", "green3") )
dev.off()



CairoPS("1a.ps")
    plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("red3", "yellow3", "cyan3"), col=c("red3", "yellow3", "cyan3") )
dev.off()
CairoPS("1b.ps")
    plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("purple3", "orange3", "green3"), col=c("purple3", "orange3", "green3") )
dev.off()



CairoSVG("1a.Cairo.svg", bg = "white")
    plot(VennDiag1, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("red3", "yellow3", "cyan3"), col=c("red3", "yellow3", "cyan3") )
dev.off()
CairoSVG("1b.Cairo.svg", bg = "white")
    plot(VennDiag2, counts = TRUE, font=1, cex=1, alpha=0.5, fill=c("purple3", "orange3", "green3"), col=c("purple3", "orange3", "green3") )
dev.off()





