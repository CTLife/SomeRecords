#########################################################################
library(methylKit)
library(genomation)




myTest_1 <- "1-CpG/3_PD10-TFC-hMII3_Rep3.CpG.bismark.cov"
myTest_2 <- "1-CpG/3_PD10-TFC-hMII4_Rep4.CpG.bismark.cov"
myCtrl_1 <- "1-CpG/6_sMII-5ge-8-27_Rep1.CpG.bismark.cov"
myCtrl_2 <- "1-CpG/8_7M2-11-18_Rep1.CpG.bismark.cov"
myOutDir <- "z"
myMinReads <- 5

myOutDir <- paste(myOutDir, "_minReads=", myMinReads, sep="")
if( ! file.exists(myOutDir) ) { dir.create(myOutDir, recursive = TRUE) }



## ----message=FALSE-------------------------------------------------------
file.list=list(myTest_1, myTest_2,  myCtrl_1, myCtrl_2 )



# read the files to a methylRawList object: myobj
myobj=methRead(file.list,
           sample.id=list("test1","test2","ctrl1","ctrl2"),
           assembly="hg38",
           treatment=c(1,1,0,0),
           context="CpG",
           pipeline = "bismarkCoverage",
           mincov = myMinReads,       ## >= n
           header = FALSE
           )


sink( file=paste(myOutDir, "1-all-the-4files.txt", sep="/") )
  print(file.list)
  print("#########################")
  print(myobj)
sink()
 
pdf( file=paste(myOutDir, "2-MethylationStats-the-4files.pdf", sep="/")  )
  getMethylationStats(myobj[[1]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  getMethylationStats(myobj[[2]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  getMethylationStats(myobj[[3]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  getMethylationStats(myobj[[4]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
dev.off()


pdf( file=paste(myOutDir, "3-CoverageStats-the-4files.pdf", sep="/")  )
  getCoverageStats(myobj[[1]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  getCoverageStats(myobj[[2]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  getCoverageStats(myobj[[3]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
  getCoverageStats(myobj[[4]], plot=TRUE, both.strands=FALSE, chunk.size=2e8)
dev.off()




## ------------------------------------------------------------------------
filtered.myobj=filterByCoverage(myobj,  lo.count=myMinReads,  lo.perc=NULL,  hi.count=NULL, hi.perc=99.99)
sink( file=paste(myOutDir, "4-all-the-4files.txt", sep="/") )
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
  #getCorrelation(meth, plot=TRUE, method="kendall",  nrow=2e8)
  #getCorrelation(meth, plot=TRUE, method="spearman", nrow=2e8)
dev.off()



pdf( file=paste(myOutDir, "6-clusterSamples.pdf", sep="/")  )
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="ward", plot=TRUE)
clusterSamples(meth, dist="correlation", method="complete", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="complete", plot=TRUE)
clusterSamples(meth, dist="correlation", method="centroid", plot=TRUE)
clusterSamples(meth, dist="euclidean",   method="centroid", plot=TRUE)
dev.off()



## ------------------------------------------------------------------------
PCASamples(meth, screeplot=TRUE)

## ------------------------------------------------------------------------
PCASamples(meth)



## ------------------------------------------------------------------------
# make some batch data frame
# this is a bogus data frame
# we don't have batch information
# for the example data
sampleAnnotation=data.frame(batch_id=c("a","a","b","b"),
                            age=c(19,34,23,40))

as=assocComp(mBase=meth,sampleAnnotation)
as

# construct a new object by removing the first pricipal component
# from percent methylation value matrix
newObj=removeComp(meth,comp=1)

## ------------------------------------------------------------------------
mat=percMethylation(meth)

# do some changes in the matrix
# this is just a toy example
# ideally you want to correct the matrix
# for batch effects
mat[mat==100]=80
 
# reconstruct the methylBase from the corrected matrix
newobj=reconstruct(mat,meth)





## ------------------------------------------------------------------------
myDiff=calculateDiffMeth(meth)

## ------------------------------------------------------------------------
# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
#
# get hypo methylated bases
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
#
#
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)

## ------------------------------------------------------------------------
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)

## ----eval=FALSE----------------------------------------------------------
#  
#  sim.methylBase1<-dataSim(replicates=6,sites=1000,
#                           treatment=c(rep(1,3),rep(0,3)),
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
#                          )
#  
#  my.diffMeth<-calculateDiffMeth(sim.methylBase1[1:,],
#                                  overdispersion="MN",test="Chisq",mc.cores=1)

## ----eval=FALSE----------------------------------------------------------
#  
#  covariates=data.frame(age=c(30,80,34,30,80,40))
#  sim.methylBase<-dataSim(replicates=6,sites=1000,
#                          treatment=c(rep(1,3),rep(0,3)),
#                          covariates=covariates,
#                          sample.ids=c(paste0("test",1:3),paste0("ctrl",1:3))
#                          )
#  my.diffMeth3<-calculateDiffMeth(sim.methylBase,
#                                  covariates=covariates,
#                                  overdispersion="MN",test="Chisq",mc.cores=1)

## ---- eval=FALSE---------------------------------------------------------
#  myDiff=calculateDiffMeth(meth,num.cores=2)

## ------------------------------------------------------------------------
library(genomation)

# read the gene BED file
gene.obj=readTranscriptFeatures(system.file("extdata", "refseq.hg18.bed.txt", 
                                           package = "methylKit"))
#
# annotate differentially methylated CpGs with 
# promoter/exon/intron using annotation data
#
annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

## ------------------------------------------------------------------------
# read the shores and flanking regions and name the flanks as shores 
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", "cpgi.hg18.bed.txt", 
                                        package = "methylKit"),
                           feature.flank.name=c("CpGi","shores"))
#
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p,"GRanges"),
                                    cpg.obj$CpGi,cpg.obj$shores,
                         feature.name="CpGi",flank.name="shores")

## ------------------------------------------------------------------------
promoters=regionCounts(myobj,gene.obj$promoters)

head(promoters[[1]])

## ------------------------------------------------------------------------
diffAnn=annotateWithGeneParts(as(myDiff25p,"GRanges"),gene.obj)

# target.row is the row number in myDiff25p
head(getAssociationWithTSS(diffAnn))

## ------------------------------------------------------------------------
getTargetAnnotationStats(diffAnn,percentage=TRUE,precedence=TRUE)

## ------------------------------------------------------------------------
plotTargetAnnotation(diffAnn,precedence=TRUE,
    main="differential methylation annotation")

## ------------------------------------------------------------------------
plotTargetAnnotation(diffCpGann,col=c("green","gray","white"),
       main="differential methylation annotation")

## ------------------------------------------------------------------------
getFeatsWithTargetsStats(diffAnn,percentage=TRUE)

## ------------------------------------------------------------------------
class(meth)
as(meth,"GRanges")
class(myDiff)
as(myDiff,"GRanges")

## ------------------------------------------------------------------------
class(myobjDB[[1]])

## ----eval=FALSE----------------------------------------------------------
#  as(myobjDB[[1]],"methylRaw")

## ----eval=FALSE----------------------------------------------------------
#  data(methylKit)
#  
#  objDB=makeMethylDB(methylBase.obj,"exMethylDB")
#  

## ------------------------------------------------------------------------
select(meth,1:5) # get first 10 rows of a methylBase object
myDiff[21:25,] # get 5 rows of a methylDiff object

## ----message=FALSE,warning=FALSE,eval=FALSE------------------------------
#  library(GenomicRanges)
#  my.win=GRanges(seqnames="chr21",
#                 ranges=IRanges(start=seq(from=9764513,by=10000,length.out=20),width=5000) )
#  
#  # selects the records that lie inside the regions
#  selectByOverlap(myobj[[1]],my.win)

## ----eval=FALSE----------------------------------------------------------
#  # creates a new methylRawList object
#  myobj2=reorganize(myobj,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )
#  # creates a new methylBase object
#  meth2 =reorganize(meth,sample.ids=c("test1","ctrl2"),treatment=c(1,0) )

## ---- eval=FALSE---------------------------------------------------------
#  # creates a matrix containing percent methylation values
#  perc.meth=percMethylation(meth)

## ---- eval=FALSE---------------------------------------------------------
#   download.file("https://dl.dropboxusercontent.com/u/1373164/H1.chr21.chr22.rds",
#                 destfile="H1.chr21.chr22.rds",method="curl")
#  
#   mbw=readRDS("H1.chr21.chr22.rds")
#  
#   # it finds the optimal number of componets as 6
#   res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10)
#  
#   # however the BIC stabilizes after 4, we can also try 4 componets
#   res=methSeg(mbw,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:4)
#  
#   # get segments to BED file
#   methSeg2bed(res,filename="H1.chr21.chr22.trial.seg.bed")
#  
#  

## ------------------------------------------------------------------------
sessionInfo() 

## ----eval=TRUE,echo=FALSE------------------------------------------------
# tidy up                  
rm(myobjDB)              
unlink(list.files(pattern = "methylDB",full.names = TRUE),recursive = TRUE)
