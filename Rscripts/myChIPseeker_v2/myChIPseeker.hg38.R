## example: Rscript  myChIPseeker.hg38.R  4-overlap-pups-rename     ART-E18-H3K4me1/1_ART-E18C-H3K4me1.bed    4-overlap-pups-annotation   

args <- commandArgs(TRUE)
print("args: ")
print(args[1])         
print(args[2])
print(args[3])         
print("#############")

inputDir   = args[1];     ## bed file in this folder (input dir)
inputBED   = args[2];     ## bed file  (input file)
outPath    = args[3];     ## output file path

# inputDir   = "4-overlap-pups-rename"
# inputBED   = "ART-E18-H3K4me1/2_ART-E18D-H3K4me1.bed"
# outPath    = "4-overlap-pups-annotation" 

inputFILE  = paste(inputDir, inputBED , sep="/")  
outPath2   = paste(outPath,  inputBED,  sep="/")      ##output file path     
if( ! file.exists(outPath)  ) { dir.create(outPath,  recursive = TRUE)  }
if( ! file.exists(outPath2) ) { dir.create(outPath2, recursive = TRUE)  }





###########################################################################
## ----style, echo=FALSE, results='asis', message=TRUE--------------------
BiocStyle::markdown()
knitr::opts_chunk$set(tidy=FALSE,   warning=TRUE,  message=TRUE)

## ----echo=FALSE, results='hide', message=TRUE---------------------------
library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ChIPseeker)

## -------------------------- loading packages ---------------------------- 
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
require(clusterProfiler)

## ------------------------------------------------------------------------
myCurrentFile <- inputFILE
print(myCurrentFile) 
peak <- readPeakFile( myCurrentFile )  ## table header must be removed
peak

## ----fig.height=20, fig.width=10------------------------------------------
pdf(file=paste(outPath2, "1_ChIP_Peaks_over_Chromosomes.pdf", sep="/"),  width=20, height=20)
covplot( peak, weightCol="V5" , lower = 0.01  )
dev.off()

pdf(file=paste(outPath2, "1_ChIP_Peaks_over_Chromosomes.chromXY.pdf", sep="/"),  width=20, height=10)
covplot( peak, weightCol="V5" , chrs=c("chrX", "chrY") , lower = 0.01    )
dev.off()


pdf(file=paste(outPath2, "1_ChIP_Peaks_over_Chromosomes.chromLength.pdf", sep="/"),  width=20, height=20)
covplot( peak, weightCol="V5" , lower = -1  )
dev.off()


## ------------------------------------------------------------------------
promoter  <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000)
promoter
tagMatrix <- getTagMatrix(peak, windows=promoter)
dim(tagMatrix)
tagMatrix[1:10, 1:10]

write.table(x=tagMatrix, file =paste(outPath2, "2_peakMatrix_TSS-10kb.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

## ----fig.cap="Heatmap of ChIP peaks binding to TSS regions", fig.align="center", fig.height=12, fig.width=4----
pdf(file = paste(outPath2, "2_Heatmap_of_ChIP_peaks_binding_to_TSS_regions.10kb.pdf", sep="/"),  width=10, height=20)
tagHeatmap(tagMatrix, xlim=c(-10000, 10000), color="red")
dev.off()

## ----eval=TRUE, fig.cap="Average Profile of ChIP peaks binding to TSS region", fig.align="center", fig.height=4, fig.width=7----
pdf(file=paste(outPath2, "3_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.10kb.pdf", sep="/"),  width=12, height=10)
plotAvgProf(tagMatrix, xlim=c(-10000, 10000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" )
dev.off()

pdf(file=paste(outPath2, "3_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.10kb.conf.pdf", sep="/"),  width=12, height=10)
plotAvgProf(tagMatrix, xlim=c(-10000, 10000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" , conf=0.95, resample=1000)
dev.off()


#####################################
promoter2  <- getPromoters(TxDb=txdb, upstream=4000, downstream=4000)
promoter2
tagMatrix2 <- getTagMatrix(peak, windows=promoter2)
dim(tagMatrix2)
tagMatrix2[1:10, 1:10]

write.table(x=tagMatrix2, file =paste(outPath2, "2_peakMatrix_TSS-4kb.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

## ----fig.cap="Heatmap of ChIP peaks binding to TSS regions", fig.align="center", fig.height=12, fig.width=4----
pdf(file = paste(outPath2, "2_Heatmap_of_ChIP_peaks_binding_to_TSS_regions.4kb.pdf", sep="/"),  width=10, height=20)
tagHeatmap(tagMatrix2, xlim=c(-4000, 4000), color="red")
dev.off()

## ----eval=TRUE, fig.cap="Average Profile of ChIP peaks binding to TSS region", fig.align="center", fig.height=4, fig.width=7----
pdf(file=paste(outPath2, "3_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.4kb.pdf", sep="/"),  width=12, height=10)
plotAvgProf(tagMatrix2, xlim=c(-4000, 4000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" )
dev.off()

pdf(file=paste(outPath2, "3_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.4kb.conf.pdf", sep="/"),  width=12, height=10)
plotAvgProf(tagMatrix2, xlim=c(-4000, 4000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" , conf=0.95, resample=1000)
dev.off()



## ------------------------------------------------------------------------
peakAnno <- annotatePeak( myCurrentFile, tssRegion=c(-3000, 3000),  TxDb=txdb, annoDb="org.Hs.eg.db"  )   
write.table(x=peakAnno, file = paste(outPath2, "4_peakAnno.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

## ----fig.cap="Genomic Annotation by pieplot", fig.align="center", fig.height=6, fig.width=8----
pdf(file=paste(outPath2, "4_Genomic_Annotation_by_pieplot.pdf", sep="/"),  width=7, height=5)
plotAnnoPie(peakAnno)
dev.off()

svg(filename=paste(outPath2, "4_Genomic_Annotation_by_pieplot.svg", sep="/"),  width=7, height=5)
plotAnnoPie(peakAnno)
dev.off()

## ----fig.cap="Genomic Annotation by barplot", fig.align="center", fig.height=4, fig.width=10----
pdf(file=paste(outPath2, "5_Genomic_Annotation_by_barplot.pdf", sep="/"),  width=7, height=3)
plotAnnoBar(peakAnno)
dev.off()

svg(filename=paste(outPath2, "5_Genomic_Annotation_by_barplot.svg", sep="/"),  width=7, height=3)
plotAnnoBar(peakAnno)
dev.off()

## ----fig.cap="Genomic Annotation by vennpie", fig.align="center", fig.height=8, fig.width=11----
pdf(file=paste(outPath2, "6_Genomic_Annotation_by_vennpie.pdf", sep="/"),  width=7, height=5)
vennpie(peakAnno)
dev.off()

svg(filename=paste(outPath2, "6_Genomic_Annotation_by_vennpie.svg", sep="/"),  width=7, height=5)
vennpie(peakAnno)
dev.off()

## ----eval=F, fig.cap="Genomic Annotation by upsetplot", fig.align="center", fig.height=8, fig.width=12----
pdf(file=paste(outPath2, "7_Genomic_Annotation_by_upsetplot.pdf", sep="/"),  width=8, height=5)
upsetplot(peakAnno)
dev.off()

svg(filename=paste(outPath2, "7_Genomic_Annotation_by_upsetplot.svg", sep="/"),  width=8, height=5)
upsetplot(peakAnno)
dev.off()

## ----eval=F, fig.cap="Genomic Annotation by upsetplot", fig.align="center", fig.height=8, fig.width=12----
pdf(file=paste(outPath2, "8_Genomic_Annotation_by_upsetplot_vennpie.pdf", sep="/"),  width=8, height=5)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

svg(filename=paste(outPath2, "8_Genomic_Annotation_by_upsetplot_vennpie.svg", sep="/"),  width=8, height=5)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

## ----fig.cap="Distribution of Binding Sites", fig.align="center", fig.height=2, fig.width=6----
pdf(file=paste(outPath2, "9_Distribution_of_peaks_to_TSS.pdf", sep="/"),  width=8, height=2)
plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS")
dev.off()

svg(filename=paste(outPath2, "9_Distribution_of_peaks_to_TSS.svg", sep="/"),  width=800, height=200)
plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS")
dev.off()
################################################################################################################################












library("ggplot2") 
 
MyTheme_1 <- function(textSize1=14, hjust1=NULL, vjust1=NULL,  angle1=NULL) {    # "hjust=1, vjust=1, angle=30" for some boxplots.
  theme(  
    line  = element_line(colour="black",  size=1.0,   linetype=1,      lineend=NULL),                                                                                        ## all line elements.          局部优先总体,下面3个也是,只对非局部设置有效.   所有线属性.
    rect  = element_rect(colour="black",  size=1.0,   linetype=1,      fill="transparent" ),                                                                                 ## all rectangluar elements.    hjust=1: 靠右对齐.   所有矩形区域属性.
    text  = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all text elements.           "serif" for a serif font. 所有文本相关属性.
    title = element_text(family="serif",  face="plain",  colour="black",  size=textSize1, hjust=0.5, vjust=0.5,   angle=0, lineheight=1.0,  margin = NULL, debug = NULL),    ## all title elements: plot, axes, legends.    hjust:水平对齐的方向.  所有标题属性.
    ## aspect.ratio = 1,   ##高宽比
    
    axis.title    = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## label of axes (element_text; inherits from text).  horizontal: 水平的, 水平线 
    axis.title.x  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis label (element_text; inherits from axis.title)
    axis.title.y  = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=90,      lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis label (element_text; inherits from axis.title)
    axis.text     = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## tick labels along axes (element_text; inherits from text). 坐标轴刻度的标签的属性.                                                         
    axis.text.x   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=hjust1, vjust=vjust1, angle=angle1,  lineheight=1.0,  margin = NULL, debug = NULL),       ## x axis tick labels (element_text; inherits from axis.text)
    axis.text.y   = element_text(family="serif", face="plain", colour="black", size=textSize1,    hjust=0.5,    vjust=0.5,    angle=0,       lineheight=1.0,  margin = NULL, debug = NULL),       ## y axis tick labels (element_text; inherits from axis.text)
    
    axis.ticks        = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## tick marks along axes (element_line; inherits from line). 坐标轴刻度线.
    axis.ticks.x      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## x axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.y      = element_line(colour="black", size=0.5, linetype=1, lineend=NULL),          ## y axis tick marks (element_line; inherits from axis.ticks)
    axis.ticks.length = grid::unit(2.0,   "mm",   data=NULL),                                      ## length of tick marks (unit), ‘"mm"’ Millimetres.  10 mm = 1 cm.  刻度线长度
    axis.line         = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## lines along axes (element_line; inherits from line). 坐标轴线
    axis.line.x       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL), 	 ## line along x axis (element_line; inherits from axis.line)
    axis.line.y       = element_line(colour="transparent", size=0.3, linetype=1, lineend=NULL),	   ## line along y axis (element_line; inherits from axis.line)
    
    legend.background    = element_rect(colour="transparent", size=1, linetype=1, fill="transparent" ), 	      ## background of legend (element_rect; inherits from rect)
    legend.spacing       = grid::unit(1, "mm", data=NULL), 	                                                    ## extra space added around legend (unit). linetype=1指的是矩形边框的类型.
    legend.key           = element_rect(colour="transparent", size=2, linetype=1, fill="transparent" ), 	      ## background underneath legend keys. 图例符号. size=1指的是矩形边框的大小.
    legend.key.size      = grid::unit(6,   "mm", data=NULL) , 	                                                ## size of legend keys   (unit; inherits from legend.key.size)
    legend.key.height    = grid::unit(6.5, "mm", data=NULL) , 	                                                ## key background height (unit; inherits from legend.key.size)
    legend.key.width     = grid::unit(8,   "mm", data=NULL) ,                                                   ## key background width  (unit; inherits from legend.key.size)
    legend.text          = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	##legend item labels. 图例文字标签.
    legend.text.align    = 0, 	                    ## alignment of legend labels (number from 0 (left) to 1 (right))
    legend.title         = element_blank(),   	    ## title of legend (element_text; inherits from title)
    legend.title.align   = 0, 	                    ## alignment of legend title (number from 0 (left) to 1 (right))
    legend.position      = "right", 	              ## the position of legends. ("left", "right", "bottom", "top", or two-element numeric vector)
    legend.direction     = "vertical",        	    ## layout of items in legends  ("horizontal" or "vertical")   图例排列方向
    legend.justification = "center",      	        ## anchor point for positioning legend inside plot ("center" or two-element numeric vector)  图例居中方式
    legend.box           = NULL, 	                  ## arrangement of multiple legends ("horizontal" or "vertical")  多图例的排列方式
    legend.box.just      = NULL, 	                  ## justification of each legend within the overall bounding box, when there are multiple legends ("top", "bottom", "left", or "right")  多图例的居中方式
    
    panel.background   = element_rect(colour="transparent", size=0.0, linetype=1, fill="transparent" ),   	## background of plotting area, drawn underneath plot (element_rect; inherits from rect)
    panel.border       = element_rect(colour="black", size=0.5, linetype=1, fill=NA ), 	                    ## border around plotting area, drawn on top of plot so that it covers tick marks and grid lines. This should be used with fill=NA (element_rect; inherits from rect)
    panel.spacing      = grid::unit(1, "mm", data=NULL) , 	                                                ## margin around facet panels (unit)  分面绘图区之间的边距
    panel.spacing.x    = grid::unit(1, "mm", data=NULL) ,
    panel.spacing.y    = grid::unit(1, "mm", data=NULL) ,
    panel.grid         = element_blank(), 	                                                                ## grid lines (element_line; inherits from line)  绘图区网格线
    panel.grid.major   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## major grid lines (element_line; inherits from panel.grid)  主网格线
    panel.grid.minor   = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## minor grid lines (element_line; inherits from panel.grid)  次网格线
    panel.grid.major.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) , 	    ## vertical major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.major.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal major grid lines (element_line; inherits from panel.grid.major)
    panel.grid.minor.x = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## vertical minor grid lines (element_line; inherits from panel.grid.minor)
    panel.grid.minor.y = element_line(colour="transparent", size=NULL, linetype=NULL, lineend=NULL) ,     	## horizontal minor grid lines (element_line; inherits from panel.grid.minor)
    
    plot.background	= element_rect(colour="transparent", size=NULL, linetype=NULL, fill="transparent" ),                                                ## background of the entire plot (element_rect; inherits from rect)  整个图形的背景
    plot.title      = element_text(family="serif", face=NULL, colour="black", size=textSize1, hjust=0.5, vjust=0.5,   angle=NULL, lineheight=NULL),     ## plot title (text appearance) (element_text; inherits from title)  图形标题
    plot.margin     = grid::unit(c(5, 5, 5, 5), "mm", data=NULL), 	                                                                                    ## margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
    
    strip.background = element_rect(colour=NULL,    size=NULL, linetype=NULL, fill=NULL ), 	                                                      ## background of facet labels (element_rect; inherits from rect)  分面标签背景
    strip.text       = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels (element_text; inherits from text)
    strip.text.x     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL), 	    ## facet labels along horizontal direction (element_text; inherits from strip.text)
    strip.text.y     = element_text(family="serif", face=NULL, colour=NULL, size=NULL, hjust=NULL, vjust=NULL, angle=NULL, lineheight=NULL)   	  ## facet labels along vertical direction (element_text; inherits from strip.text) 
  ) 
} 

MySaveGgplot2_1 <- function(path1, fileName1,  height1, width1) {
  SVG1 <- paste(path1,  "/",  "SVG",  sep = "",  collapse = NULL)
  PNG1 <- paste(path1,  "/",  "PNG",  sep = "",  collapse = NULL)
  PDF1 <- paste(path1,  "/",  "PDF",  sep = "",  collapse = NULL)
  EPS1 <- paste(path1,  "/",  "EPS",  sep = "",  collapse = NULL)
  if( ! file.exists(SVG1) ) { dir.create(SVG1) }
  if( ! file.exists(PNG1) ) { dir.create(PNG1) }
  if( ! file.exists(PDF1) ) { dir.create(PDF1) }
  if( ! file.exists(EPS1) ) { dir.create(EPS1) }
  ggsave( filename = paste(SVG1,  "/",  fileName1,  ".svg",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PNG1,  "/",  fileName1,  ".png",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(PDF1,  "/",  fileName1,  ".pdf",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200 )
  ggsave( filename = paste(EPS1,  "/",  fileName1,  ".eps",  sep="",  collapse=NULL),     height=height1,    width=width1,      dpi = 1200,   device=cairo_ps)         
}






## ------------------------------------------------------------------------
## DOSE for Disease Ontology. 
## ReactomePA for reactome pathway. 
## clusterProfiler for Gene Ontology and KEGG enrichment analysis.

library(DOSE)
library(ReactomePA)
library(clusterProfiler)  ## enrichGO for GO Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment GO categories after FDR control.

bp1 <- enrichGO(as.data.frame(peakAnno)$geneId,  OrgDb='org.Hs.eg.db',  ont="BP",  readable=TRUE)
dim(bp1)
bp1[1:10,1:7]

write.table(x= bp1 , file = paste(outPath2, "10_BP_enrichment.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n",  na = "NA",   dec = ".",   row.names = TRUE,
            col.names = TRUE,   qmethod = c("escape", "double"),
            fileEncoding = "")


pdf(file=paste(outPath2, "10A_BP.pdf", sep="/"),  width=8, height=5)
dotplot(bp1, showCategory = 20)
dev.off()

png(filename=paste(outPath2, "10A_BP.png", sep="/"),  width=800, height=500)
dotplot(bp1, showCategory = 20)
dev.off()

svg(filename=paste(outPath2, "10A_BP.svg", sep="/"),  width=8, height=5)
dotplot(bp1, showCategory = 20)
dev.off()

##################################
bp2     <- as.data.frame( bp1 )   
dim(bp2)
bp2[1:20, -8]
bp2$ID

bp2$GeneRatio <- sapply(bp1$GeneRatio, function(x) eval(parse(text=x))) * 100
bp2$pvalue    <- -log10(bp2$pvalue)
bp2$p.adjust  <- -log10(bp2$p.adjust)
bp2$qvalue    <- -log10(bp2$qvalue)

bp3 <- bp2[1:20, -8]
bp3
bp3 <- transform(bp3, ID= rev( factor(ID, levels=unique(ID))) )
bp3
myMidValue <- ( max(bp3$GeneRatio) + min(bp3$GeneRatio)  )/2
FigureTemp1 <- ggplot( data = bp3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
               geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
               scale_y_discrete( breaks=bp3$ID, labels=rev(bp3$ID) ) + 
               xlab("-log10(ajusted p-value)") +   ylab("GO terms (BP)") +   ggtitle("GO terms (BP)") +  theme_bw()
MySaveGgplot2_1( path1=outPath2, fileName1="10B_BP",  height1=6,  width1=5)





#####################################################################
MF1 <- enrichGO(as.data.frame(peakAnno)$geneId,  OrgDb='org.Hs.eg.db',  ont="MF",  readable=TRUE)
dim(MF1)
MF1[1:10,1:7]

write.table(x= MF1 , file = paste(outPath2, "10_MF_enrichment.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n",  na = "NA",   dec = ".",   row.names = TRUE,
            col.names = TRUE,   qmethod = c("escape", "double"),
            fileEncoding = "")


pdf(file=paste(outPath2, "10A_MF.pdf", sep="/"),  width=8, height=5)
dotplot(MF1, showCategory = 20)
dev.off()

png(filename=paste(outPath2, "10A_MF.png", sep="/"),  width=800, height=500)
dotplot(MF1, showCategory = 20)
dev.off()

svg(filename=paste(outPath2, "10A_MF.svg", sep="/"),  width=8, height=5)
dotplot(MF1, showCategory = 20)
dev.off()

##################
MF2     <- as.data.frame( MF1 )   
dim(MF2)
MF2[1:10, -8]
MF2$ID

MF2$GeneRatio <- sapply(MF1$GeneRatio, function(x) eval(parse(text=x))) * 100
MF2$pvalue    <- -log10(MF2$pvalue)
MF2$p.adjust  <- -log10(MF2$p.adjust)
MF2$qvalue    <- -log10(MF2$qvalue)

MF3 <- MF2[1:20, -8]
MF3 <- transform(MF3, ID= rev( factor(ID, levels=unique(ID))) )
myMidValue <- ( max(MF3$GeneRatio) + min(MF3$GeneRatio)  )/2
FigureTemp1 <- ggplot( data = MF3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
    geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
    scale_y_discrete( breaks=MF3$ID, labels=rev(MF3$ID) ) + 
    xlab("-log10(ajusted p-value)") +   ylab("GO terms (MF)") +   ggtitle("GO terms (MF)") + 
    theme_bw()
MySaveGgplot2_1( path1=outPath2, fileName1="10B_MF",  height1=6,  width1=5)





#####################################################################################################
CC1 <- enrichGO(as.data.frame(peakAnno)$geneId,  OrgDb='org.Hs.eg.db',  ont="CC",  readable=TRUE)
dim(CC1)
CC1[1:10,1:7]

write.table(x= CC1 , file = paste(outPath2, "10_CC_enrichment.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n",  na = "NA",   dec = ".",   row.names = TRUE,
            col.names = TRUE,   qmethod = c("escape", "double"),
            fileEncoding = "")


pdf(file=paste(outPath2, "10A_CC.pdf", sep="/"),  width=8, height=5)
dotplot(CC1, showCategory = 20)
dev.off()

png(filename=paste(outPath2, "10A_CC.png", sep="/"),  width=800, height=500)
dotplot(CC1, showCategory = 20)
dev.off()

svg(filename=paste(outPath2, "10A_CC.svg", sep="/"),  width=8, height=5)
dotplot(CC1, showCategory = 20)
dev.off()

#####################
CC2     <- as.data.frame( CC1 )   
dim(CC2)
CC2[1:10, -8]
CC2$ID

CC2$GeneRatio <- sapply(CC1$GeneRatio, function(x) eval(parse(text=x))) * 100
CC2$pvalue    <- -log10(CC2$pvalue)
CC2$p.adjust  <- -log10(CC2$p.adjust)
CC2$qvalue    <- -log10(CC2$qvalue)

CC3 <- CC2[1:20, -8]
CC3 <- transform(CC3, ID= rev( factor(ID, levels=unique(ID))) )
myMidValue <- ( max(CC3$GeneRatio) + min(CC3$GeneRatio)  )/2
FigureTemp1 <- ggplot( data = CC3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
  geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
  scale_y_discrete( breaks=CC3$ID, labels=rev(CC3$ID) ) + 
  xlab("-log10(ajusted p-value)") +   ylab("GO terms (CC)") +   ggtitle("GO terms (CC)") + 
  theme_bw()
MySaveGgplot2_1( path1=outPath2, fileName1="10B_CC",  height1=6,  width1=5)





################################################################################################################################
library(ReactomePA)
length( as.data.frame(peakAnno)$geneId )
pathway1 <- enrichPathway(gene=as.data.frame(peakAnno)$geneId, organism = "human")   ## enrichPathway {ReactomePA}
dim(pathway1)
head(pathway1, 2)

pdf(file=paste(outPath2, "11A_pathway1.pdf", sep="/"),  width=8, height=5)
dotplot(pathway1, showCategory = 20)
dev.off()

png(filename=paste(outPath2, "11A_pathway1.png", sep="/"),  width=800, height=500)
dotplot(pathway1, showCategory = 20)
dev.off()

svg(filename=paste(outPath2, "11A_pathway1.svg", sep="/"),  width=8, height=5)
dotplot(pathway1, showCategory = 20)
dev.off()


write.table(x= pathway1 , file = paste(outPath2, "11_pathway1.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n",  na = "NA",   dec = ".",   row.names = TRUE,
            col.names = TRUE,   qmethod = c("escape", "double"),
            fileEncoding = "")


pdf(file=paste(outPath2, "11A_pathway1.pdf", sep="/"),  width=8, height=5)
dotplot(pathway1, showCategory = 20)
dev.off()

png(filename=paste(outPath2, "11A_pathway1.png", sep="/"),  width=800, height=500)
dotplot(pathway1, showCategory = 20)
dev.off()

svg(filename=paste(outPath2, "11A_pathway1.svg", sep="/"),  width=8, height=5)
dotplot(pathway1, showCategory = 20)
dev.off()

pathway12     <- as.data.frame( pathway1 )   
dim(pathway12)
pathway12[1:10, -8]
pathway12$ID

pathway12$GeneRatio <- sapply(pathway12$GeneRatio, function(x) eval(parse(text=x))) * 100
pathway12$pvalue    <- -log10(pathway12$pvalue)
pathway12$p.adjust  <- -log10(pathway12$p.adjust)
pathway12$qvalue    <- -log10(pathway12$qvalue)

pathway13 <- pathway12[1:20, ]
pathway13 <- transform(pathway13, ID= rev( factor(ID, levels=unique(ID))) )
myMidValue <- ( max(pathway13$GeneRatio) + min(pathway13$GeneRatio)  )/2
FigureTemp1 <- ggplot( data = pathway13, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
  geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
  xlab("-log10(ajusted p-value)") +   ylab("pathway") +   ggtitle("pathway") +  scale_y_discrete( breaks=pathway13$ID, labels=rev(pathway13$ID) ) + 
  theme_bw()
MySaveGgplot2_1(   path1=outPath2, fileName1="11B_pathway1",  height1=6,  width1=5)



###########################################################################################
gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
length(gene)
gene[1:10]
pathway2 <- enrichPathway(gene, organism = "human")
dim(pathway2)
head(pathway2, 2)
write.table(x= pathway2 , file = paste(outPath2, "12_pathway2.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t",  eol = "\n",  na = "NA",   dec = ".",  
            row.names = TRUE,  col.names = TRUE,   qmethod = c("escape", "double"), fileEncoding = "")

pdf(file=paste(outPath2, "12A_pathway2.pdf", sep="/"),  width=8, height=5)
dotplot(pathway2, showCategory = 20)
dev.off()

svg(filename=paste(outPath2, "12A_pathway2.svg", sep="/"),  width=8, height=5)
dotplot(pathway2, showCategory = 20)
dev.off()


pathway22     <- as.data.frame( pathway2 )   
dim(pathway22)
pathway22[1:10, -8]
pathway22$ID

pathway22$GeneRatio <- sapply(pathway22$GeneRatio, function(x) eval(parse(text=x))) * 100
pathway22$pvalue    <- -log10(pathway22$pvalue)
pathway22$p.adjust  <- -log10(pathway22$p.adjust)
pathway22$qvalue    <- -log10(pathway22$qvalue)

pathway23 <- pathway22[1:20, ]
pathway23 <- transform(pathway23, ID= rev( factor(ID, levels=unique(ID))) )
myMidValue <- ( max(pathway23$GeneRatio) + min(pathway23$GeneRatio)  )/2
FigureTemp1 <- ggplot( data = pathway23, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
  geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
  xlab("-log10(ajusted p-value)") +   ylab("pathway") +   ggtitle("pathway") + scale_y_discrete( breaks=pathway23$ID, labels=rev(pathway23$ID) ) + 
  theme_bw()
MySaveGgplot2_1(  path1=outPath2, fileName1="12B_pathway2", height1=6, width1=5)




################################################################################
library(stringi)

all_nearest_genes   <- as.data.frame(peakAnno)$geneId
all_nearest_genes2 <- do.call(paste, c(as.list(all_nearest_genes), sep=","))  
length(all_nearest_genes)
length(all_nearest_genes2)
nearest_genes_10kb  <- all_nearest_genes[ abs( as.data.frame(peakAnno)$distanceToTSS) < 10000 ]
nearest_genes_10kb_2 <- do.call(paste, c(as.list(nearest_genes_10kb), sep=","))  
length(nearest_genes_10kb)
length(nearest_genes_10kb_2)

shell_cmd1 <- paste("Rscript   clusterProfiler.oneList.hg38.R  ",   all_nearest_genes2 ,   "    ",  outPath2,   "/all_nearest_genes" ,  sep="") 
shell_cmd2 <- paste("Rscript   clusterProfiler.oneList.hg38.R  ",   nearest_genes_10kb_2,  "    ",  outPath2,   "/nearest_genes_10kb",  sep="") 
system( shell_cmd1 )
system( shell_cmd2 )



 









