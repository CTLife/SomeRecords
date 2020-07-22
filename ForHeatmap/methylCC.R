
library(methylCC)
library(ggplot2)
library(bsseq)
library(methylSig)
library(stats)


 
  
files1 = read.table(file="AllPups.5groups.txt", header = TRUE, sep = "\t")

files1_name = files1[,1]
files1_full = paste(files1[,2], files1[,1], sep="/" )
files1_group = files1[,3]

matrix1 = read.table(file="BED/final.bed", header = TRUE, sep = "\t")
head(matrix1)

GRanges1 = makeGRangesFromDataFrame(df = matrix1,
                         keep.extra.columns=FALSE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)





BSseq_object = bsseq::read.bismark(files = files1_full,
             loci = GRanges1,
             colData = DataFrame(row.names = files1_name, type = files1_group ),
             rmZeroCov = FALSE,
             strandCollapse = TRUE,
             nThread = 4L  )


pData(BSseq_object)
rownames( pData(BSseq_object) )
colnames( pData(BSseq_object) )

unique( files1_group )
length( files1_group[files1_group=="1NC"] )
length( files1_group[files1_group=="2IVFfresh"] )
length( files1_group[files1_group=="3ICSIfresh"] )
length( files1_group[files1_group=="4IVFfrozen"] )
length( files1_group[files1_group=="5ICSIfrozen"] )

BSseq_object2 = methylSig::filter_loci_by_group_coverage(bs=BSseq_object, group_column="type", 
                 min_samples_per_group=c('1NC'=34,  '2IVFfresh'=26, '3ICSIfresh'=30, '4IVFfrozen'=28, '5ICSIfrozen'=19) )


getCoverage(BSseq_object2)
getMeth(BSseq_object2, type = "raw" )
getStats(BSseq_object2)
poissonGoodnessOfFit(BSseq_object2)
binomialGoodnessOfFit(BSseq_object2)


BSseq_object3 = methylSig::filter_loci_by_coverage(BSseq_object2, min_count = 6, max_count = 1000)
 
set.seed(12345)
est <- estimatecc(object = BSseq_object2,  include_cpgs = TRUE , include_dmrs = TRUE ) 
est
matrixC = cell_counts(est)

write.table( x=matrixC, file = "matrixC.txt", append = FALSE, quote = FALSE, sep = "\t",
             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
             col.names = TRUE )





##heatmap3 is the central function of the heatmap3 package. 
##Beware that this is different from "heatmap.3", of which there are numerous versions

library(heatmap3)
matrixC
matrixC_anno = data.frame( "name"=files1_name, "type"=files1_group)
 


# Make helper function to map metadata category to color
mapTypeToColor<-function(annotations){ 
  colorsVector = ifelse(annotations["type"]=="1NC",       "red", 
                 ifelse(annotations["type"]=="2IVFfresh",  "orange",
                 ifelse(annotations["type"]=="3ICSIfresh",  "yellow4",
                 ifelse(annotations["type"]=="4IVFfrozen",  "blue",
                 ifelse(annotations["type"]=="5ICSIfrozen",  "cyan", 'black'        
                  )))))
  return(colorsVector)
}

mapTypeToColor(matrixC_anno)



# Test heatmap3 with several annotation options
myHeatmap3_1 <- function(myMatrix, annotations) {    
  sampleColors = mapTypeToColor(annotations)
  
  # Assign just column annotations
  heatmap3(myMatrix, margins=c(5,8), ColSideColors=sampleColors) 
  
  # Assign column annotations and make a custom legend for them
  heatmap3(myMatrix, margins=c(5,8), ColSideColors=sampleColors, 
           legendfun=function()showLegend(legend=c("MiracleDrugA", 
                                                   "MiracleDrugB", "?"), col=c("blue", "green", "red"), cex=1.5))
  
  # Assign column annotations as a mini-graph instead of colors,
  # and use the built-in labeling for them
  ColSideAnn <- data.frame(Drug=annotations[["subject_drug"]])
  heatmap3(logCPM,ColSideAnn=ColSideAnn,
           #ColSideFun=function(x)showAnn(x),
           ColSideWidth=0.8)
}

testHeatmap3(logCPM=matrixC, annotations=gAnnotationData)


pdf("1.pdf")
sampleColors = mapTypeToColor(matrixC_anno)  
heatmap3(matrixC, margins=c(5,8), RowSideColors=sampleColors, balanceColor=T, scale="none" )
dev.off()









