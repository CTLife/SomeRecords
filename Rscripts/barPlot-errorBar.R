library(reshape)                      
library("sciplot")                   
TAC1 <- read.table("echo-data.txt", header=TRUE)                            

pdf(file= "TAC_echo.pdf",  width = 3, height = 5, pointsize = 1)                   
bargraph.CI(x.factor=group, response=value,    col=c("banding"="blue", "sham"="green4"), 
            angle=NULL, density=NULL,    lc=TRUE, uc=TRUE, legend=FALSE, ncol=1,    
            leg.lab=NULL, x.leg=NULL, y.leg=NULL, cex.leg=1,    bty="n", bg="white",     
            err.col="black", err.lty=1,        xlab = "banding vs sham", ylab = "FS (%)",   
            ylim=c(0,50), xpd=FALSE, data=TAC1, subset=NULL )                         
dev.off()                               
