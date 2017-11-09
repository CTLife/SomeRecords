#########################################################################
file.list=list(
"1-CpG/3_PD10-TFC-hMII3_Rep3.CpG.bismark.cov",
"1-CpG/3_PD10-TFC-hMII4_Rep4.CpG.bismark.cov",
"1-CpG/6_sMII-5ge-8-27_Rep1.CpG.bismark.cov",
"1-CpG/8_7M2-11-18_Rep1.CpG.bismark.cov"
)

mySampleID <- list("PD10_1", "PD10_2", "M2_1", "M2_2"  )
myTreatment <- c(1,  1, 2, 2 )


#########################################################################
library(methylKit)
library(genomation)

myOutDir <- "z"
myMinReads <- 5

myOutDir <- paste(myOutDir, "_minReads=", myMinReads, sep="")
if( ! file.exists(myOutDir) ) { dir.create(myOutDir, recursive = TRUE) }



# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=mySampleID,
           assembly="hg38",
           treatment=myTreatment,
           context="CpG",
           pipeline = "bismarkCoverage",
           mincov = myMinReads,       ## >= n
           header = FALSE
           )


sink( file=paste(myOutDir, "1-all-the-files.txt", sep="/") )
  print(file.list)
  print("#########################")
  print(myobj)
sink()
 




pdf( file=paste(myOutDir, "2-MethylationStats-the-files.pdf", sep="/")  )
  for( i in c(1:length(file.list)) ) {
    getMethylationStats(myobj[[i]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  }
dev.off()


pdf( file=paste(myOutDir, "3-CoverageStats-the-files.pdf", sep="/")  )
for( i in c(1:length(file.list)) ) {
  getCoverageStats(myobj[[i]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
}
dev.off()

 



## ------------------------------------------------------------------------
filtered.myobj=filterByCoverage(myobj,  lo.count=myMinReads,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.99)
sink( file=paste(myOutDir, "4-all-the-files.txt", sep="/") )
   print(filtered.myobj)
sink()



## Merging samples
## Setting destrand=TRUE (the default is FALSE) will merge reads on both strands of a CpG dinucleotide. 
## This provides better coverage, but only advised when looking at CpG methylation (for CpH methylation this will cause wrong results in subsequent analyses).
meth=unite(myobj, destrand=FALSE)
head(meth)


## ----eval=FALSE----------------------------------------------------------
#  # creates a methylBase object,
#  # where only CpGs covered with at least 1 sample per group will be returned
#  
#  # there were two groups defined by the treatment vector,
#  # given during the creation of myobj: treatment=c(1,1,0,0)
# meth.min=unite(myobj,  min.per.group=1L)


pdf( file=paste(myOutDir, "5-getCorrelation.pdf", sep="/")  )
  getCorrelation(meth, plot=TRUE, method="pearson",  nrow=2e8)
  getCorrelation(meth, plot=TRUE, method="kendall",  nrow=2e8)
  getCorrelation(meth, plot=TRUE, method="spearman", nrow=2e8)
dev.off()



pdf( file=paste(myOutDir, "6-clusterSamples.pdf", sep="/")  )
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="ward", plot=TRUE)
clusterSamples(meth, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()


library(ggfortify)


pdf( file=paste(myOutDir, "7-PCASamples-importance-of-components.pdf", sep="/")  )
PCASamples(meth, screeplot=TRUE)
dev.off()


pdf( file=paste(myOutDir, "8-PCASamples.pdf", sep="/")  )
summary( PCASamples(meth, adj.lim=c(0, 0)) )
dev.off()






