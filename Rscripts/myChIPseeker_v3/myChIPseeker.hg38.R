
######################################################################################################################################################
## example: Rscript  myChIPseeker.hg38.R    5-peaksOverlap-pups-rename     ART-E18-H3K4me1/1_ART-E18C-H3K4me1.bed    6-overlapPupsPeaks-annotation  

args <- commandArgs(TRUE)
print("args: ")
print(args[1])         
print(args[2])
print(args[3])         
print("#############")

inputDir   = args[1];     ## bed file in this folder (input dir)
inputBED   = args[2];     ## bed file  (input file)
outPath    = args[3];     ## output file path
# inputDir   = "5-peaksOverlap-pups-rename"
# inputBED   = "ART-vs-NC-H3K4me1/2_E29-H3K4me1-D-RT.bed"
# outPath    = "6-overlapPupsPeaks-annotation" 

inputFILE  = paste(inputDir, inputBED , sep="/")  
outPath2   = paste(outPath,  inputBED,  sep="/")      ##output file path     
if( ! file.exists(outPath)  ) { dir.create(outPath,  recursive = TRUE)  }
if( ! file.exists(outPath2) ) { dir.create(outPath2, recursive = TRUE)  }
######################################################################################################################################################





######################################################################################################################################################
## ----style, echo=FALSE, results='asis', message=TRUE--------------------
BiocStyle::markdown()
knitr::opts_chunk$set(tidy=FALSE,   warning=TRUE,  message=TRUE)

continue_on_error <- function()
{
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set
'options(error=continue_on_error())'")
}

# This is the key option
options(error=continue_on_error) 


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
######################################################################################################################################################





######################################################################################################################################################
my_folder_1 <- paste(outPath2,  "1-peaks-on-chromosomes",  sep="/") 
if( ! file.exists(my_folder_1)  ) { dir.create(my_folder_1,  recursive = TRUE)  }

## ----fig.height=20, fig.width=10------------------------------------------
pdf(file=paste(my_folder_1, "1A_ChIP_Peaks_over_Chromosomes.pdf", sep="/"),  width=20, height=20)
covplot( peak, weightCol="V5" , lower = 0.01  )
dev.off()

pdf(file=paste(my_folder_1, "1B_ChIP_Peaks_over_Chromosomes.chromXY.pdf", sep="/"),  width=20, height=10)
covplot( peak, weightCol="V5" , chrs=c("chrX", "chrY") , lower = 0.01    )
dev.off()

pdf(file=paste(my_folder_1, "1C_ChIP_Peaks_over_Chromosomes.chromLength.pdf", sep="/"),  width=20, height=20)
covplot( peak, weightCol="V5" , lower = -1  )
dev.off()
######################################################################################################################################################





######################################################################################################################################################
my_folder_2 <- paste(outPath2,  "2-peaks-on-TSSs-10kb",  sep="/") 
if( ! file.exists(my_folder_2)  ) { dir.create(my_folder_2,  recursive = TRUE)  }

## ------------------------------------------------------------------------
promoter  <- getPromoters(TxDb=txdb, upstream=10000, downstream=10000)
promoter
tagMatrix <- getTagMatrix(peak, windows=promoter)
dim(tagMatrix)

write.table(x=tagMatrix, file =paste(my_folder_2, "2A_peakMatrix_TSS-10kb.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

## ----fig.cap="Heatmap of ChIP peaks binding to TSS regions", fig.align="center", fig.height=12, fig.width=4----
pdf(file = paste(my_folder_2, "2B_Heatmap_of_ChIP_peaks_binding_to_TSS_regions.10kb.pdf", sep="/"),  width=10, height=20)
tagHeatmap(tagMatrix, xlim=c(-10000, 10000), color="red")
dev.off()

## ----eval=TRUE, fig.cap="Average Profile of ChIP peaks binding to TSS region", fig.align="center", fig.height=4, fig.width=7----
pdf(file=paste(my_folder_2, "2C_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.10kb.pdf", sep="/"),  width=12, height=10)
plotAvgProf(tagMatrix, xlim=c(-10000, 10000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" )
dev.off()

#pdf(file=paste(my_folder_2, "2D_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.10kb.conf.pdf", sep="/"),  width=12, height=10)
#plotAvgProf(tagMatrix, xlim=c(-10000, 10000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" , conf=0.95, resample=1000)
#dev.off()
######################################################################################################################################################





######################################################################################################################################################
my_folder_3 <- paste(outPath2,  "3-peaks-on-TSSs-4kb",  sep="/") 
if( ! file.exists(my_folder_3)  ) { dir.create(my_folder_3,  recursive = TRUE)  }

promoter2  <- getPromoters(TxDb=txdb, upstream=4000, downstream=4000)
promoter2
tagMatrix2 <- getTagMatrix(peak, windows=promoter2)
dim(tagMatrix2)

write.table(x=tagMatrix2, file =paste(my_folder_3, "3A_peakMatrix_TSS-4kb.txt", sep="/"), 
            append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
            row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

## ----fig.cap="Heatmap of ChIP peaks binding to TSS regions", fig.align="center", fig.height=12, fig.width=4----
pdf(file = paste(my_folder_3, "3B_Heatmap_of_ChIP_peaks_binding_to_TSS_regions.4kb.pdf", sep="/"),  width=10, height=20)
tagHeatmap(tagMatrix2, xlim=c(-4000, 4000), color="red")
dev.off()

## ----eval=TRUE, fig.cap="Average Profile of ChIP peaks binding to TSS region", fig.align="center", fig.height=4, fig.width=7----
pdf(file=paste(my_folder_3, "3C_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.4kb.pdf", sep="/"),  width=12, height=10)
plotAvgProf(tagMatrix2, xlim=c(-4000, 4000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" )
dev.off()

#pdf(file=paste(my_folder_3, "3D_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.4kb.conf.pdf", sep="/"),  width=12, height=10)
#plotAvgProf(tagMatrix2, xlim=c(-4000, 4000),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" , conf=0.95, resample=1000)
#dev.off()
######################################################################################################################################################





######################################################################################################################################################
my_folder_4 <- paste(outPath2,  "4-peaks-distribution",  sep="/") 
if( ! file.exists(my_folder_4)  ) { dir.create(my_folder_4,  recursive = TRUE)  }

## ------------------------------------------------------------------------
peakAnno <- annotatePeak( myCurrentFile, tssRegion=c(-3000, 3000),  TxDb=txdb, annoDb="org.Hs.eg.db"  )  
sink(  file = paste(my_folder_4, "4A_peakAnno_summary.txt", sep="/") )
print(peakAnno)  
sink()
write.table(x=peakAnno, file = paste(my_folder_4, "4B_peakAnno.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")

## ----fig.cap="Genomic Annotation by pieplot", fig.align="center", fig.height=6, fig.width=8----
pdf(file=paste(my_folder_4, "4C_Genomic_Annotation_by_pieplot.pdf", sep="/"),  width=7, height=5)
plotAnnoPie(peakAnno)
dev.off()

svg(filename=paste(my_folder_4, "4C_Genomic_Annotation_by_pieplot.svg", sep="/"),  width=7, height=5)
plotAnnoPie(peakAnno)
dev.off()

## ----fig.cap="Genomic Annotation by barplot", fig.align="center", fig.height=4, fig.width=10----
pdf(file=paste(my_folder_4, "4D_Genomic_Annotation_by_barplot.pdf", sep="/"),  width=7, height=3)
plotAnnoBar(peakAnno)
dev.off()

svg(filename=paste(my_folder_4, "4D_Genomic_Annotation_by_barplot.svg", sep="/"),  width=7, height=3)
plotAnnoBar(peakAnno)
dev.off()

## ----fig.cap="Genomic Annotation by vennpie", fig.align="center", fig.height=8, fig.width=11----
pdf(file=paste(my_folder_4, "4E_Genomic_Annotation_by_vennpie.pdf", sep="/"),  width=7, height=5)
vennpie(peakAnno)
dev.off()

svg(filename=paste(my_folder_4, "4E_Genomic_Annotation_by_vennpie.svg", sep="/"),  width=7, height=5)
vennpie(peakAnno)
dev.off()

## ----eval=F, fig.cap="Genomic Annotation by upsetplot", fig.align="center", fig.height=8, fig.width=12----
pdf(file=paste(my_folder_4, "4F_Genomic_Annotation_by_upsetplot.pdf", sep="/"),  width=8, height=5)
upsetplot(peakAnno)
dev.off()

svg(filename=paste(my_folder_4, "4F_Genomic_Annotation_by_upsetplot.svg", sep="/"),  width=8, height=5)
upsetplot(peakAnno)
dev.off()

## ----eval=F, fig.cap="Genomic Annotation by upsetplot", fig.align="center", fig.height=8, fig.width=12----
pdf(file=paste(my_folder_4, "4G_Genomic_Annotation_by_upsetplot_vennpie.pdf", sep="/"),  width=8, height=5)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

svg(filename=paste(my_folder_4, "4G_Genomic_Annotation_by_upsetplot_vennpie.svg", sep="/"),  width=8, height=5)
upsetplot(peakAnno, vennpie=TRUE)
dev.off()

## ----fig.cap="Distribution of Binding Sites", fig.align="center", fig.height=2, fig.width=6----
pdf(file=paste(my_folder_4, "4H_Distribution_of_peaks_to_TSS.pdf", sep="/"),  width=8, height=2)
plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS")
dev.off()

svg(filename=paste(my_folder_4, "4H_Distribution_of_peaks_to_TSS.svg", sep="/"),  width=800, height=200)
plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS")
dev.off()
######################################################################################################################################################





######################################################################################################################################################
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

MyGeneOntology_1 <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName="5_allNearestGenes") {
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS)<MyMaxDis, ]
  print(nrow(peakAnno2))
  
  bp1 <- enrichGO(peakAnno2$geneId,  OrgDb='org.Hs.eg.db',  ont="BP",  readable=TRUE)
  write.table(x= bp1 , file = paste(MyFolder, "/", MyFileName, "_1A_BP_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE,   qmethod = c("escape", "double"), fileEncoding = "")
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_1B_BP_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(bp1, showCategory = 20) )
  dev.off()
  
  png(filename=paste(MyFolder, "/", MyFileName, "_1B_BP_enrichment.png", sep=""),  width=800, height=500)
  print( dotplot(bp1, showCategory = 20) )
  dev.off()
  
  svg(filename=paste(MyFolder, "/", MyFileName, "_1B_BP_enrichment.svg", sep=""),  width=8, height=5)
  print( dotplot(bp1, showCategory = 20) )
  dev.off()
  
  ##################################
  bp2     <- as.data.frame( bp1 ) 
  if(nrow(bp2) > 1) {
    bp2$GeneRatio <- sapply(bp1$GeneRatio, function(x) eval(parse(text=x))) * 100
    bp2$pvalue    <- -log10(bp2$pvalue)
    bp2$p.adjust  <- -log10(bp2$p.adjust)
    bp2$qvalue    <- -log10(bp2$qvalue)
    
    bp3 <- bp2[1:20, -8]
    if( nrow(bp2) <20 ) {bp3 <- bp2[, -8]}
    bp3 <- transform(bp3, ID= rev( factor(ID, levels=unique(ID))) )

    myMidValue <- ( max(bp3$GeneRatio) + min(bp3$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = bp3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=bp3$ID, labels=rev(bp3$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("GO terms (BP)") +   ggtitle("GO terms (BP)") +  theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_1C_BP_enrichment", sep=""),  height1=6,  width1=5)
  }
  
  
  
  
  #####################################################################
  MF1 <- enrichGO(peakAnno2$geneId,  OrgDb='org.Hs.eg.db',  ont="MF",  readable=TRUE)
  write.table(x= MF1 , file = paste(MyFolder, "/", MyFileName, "_2A_MF_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE,  col.names = TRUE,   qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_2B_MF_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(MF1, showCategory = 20) )
  dev.off()
  
  png(filename=paste(MyFolder, "/", MyFileName, "_2B_MF_enrichment.png", sep=""),  width=800, height=500)
  print( dotplot(MF1, showCategory = 20) )
  dev.off()
  
  svg(filename=paste(MyFolder, "/", MyFileName, "_2B_MF_enrichment.svg", sep=""),  width=8, height=5)
  print( dotplot(MF1, showCategory = 20) )
  dev.off()
  
  ##################
  MF2     <-  as.data.frame( MF1 )  

  if(nrow(MF2) > 1) {
    MF2$GeneRatio <- sapply(MF1$GeneRatio, function(x) eval(parse(text=x))) * 100
    MF2$pvalue    <- -log10(MF2$pvalue)
    MF2$p.adjust  <- -log10(MF2$p.adjust)
    MF2$qvalue    <- -log10(MF2$qvalue)
    
    MF3 <- MF2[1:20, -8]
    if( nrow(MF2) <20 ) {MF3 <- MF2[, -8]}
    
    MF3 <- transform(MF3, ID= rev( factor(ID, levels=unique(ID))) )
    myMidValue <- ( max(MF3$GeneRatio) + min(MF3$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = MF3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=MF3$ID, labels=rev(MF3$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("GO terms (MF)") +   ggtitle("GO terms (MF)") + 
      theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_2C_MF_enrichment", sep=""),  height1=6,  width1=5)
  }
  
  
  
  
  #####################################################################################################
  CC1 <- enrichGO(peakAnno2$geneId,  OrgDb='org.Hs.eg.db',  ont="CC",  readable=TRUE)

  write.table(x= CC1 , file = paste(MyFolder, "/", MyFileName, "_3A_CC_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE,  col.names = TRUE,   qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_3B_CC_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(CC1, showCategory = 20) )
  dev.off()
  
  png(filename=paste(MyFolder, "/", MyFileName, "_3B_CC_enrichment.png", sep=""),  width=800, height=500)
  print( dotplot(CC1, showCategory = 20) )
  dev.off()
  
  svg(filename=paste(MyFolder, "/", MyFileName, "_3B_CC_enrichment.svg", sep=""),  width=8, height=5)
  print( dotplot(CC1, showCategory = 20) )
  dev.off()
  
  #####################
  CC2     <- as.data.frame( CC1  )

  if(nrow(CC2) > 1) {
    CC2$GeneRatio <- sapply(CC1$GeneRatio, function(x) eval(parse(text=x))) * 100
    CC2$pvalue    <- -log10(CC2$pvalue)
    CC2$p.adjust  <- -log10(CC2$p.adjust)
    CC2$qvalue    <- -log10(CC2$qvalue)
    
    CC3 <- CC2[1:20, -8]
    if( nrow(CC2) <20 ) {CC3 <- CC2[, -8]}
    
    CC3 <- transform(CC3, ID= rev( factor(ID, levels=unique(ID))) )
    myMidValue <- ( max(CC3$GeneRatio) + min(CC3$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = CC3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=CC3$ID, labels=rev(CC3$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("GO terms (CC)") +   ggtitle("GO terms (CC)") + 
      theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste( MyFileName, "_3C_CC_enrichment", sep=""),  height1=6,  width1=5)
  }
  
  
  
  

  
  
  
  #####################################################################################################
  Reactome1 <- enrichPathway(gene=peakAnno2$geneId, organism = "human")   ## enrichPathway {ReactomePA}
  write.table(x= Reactome1 , file = paste(MyFolder, "/", MyFileName, "_4A_Reactome_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE,  col.names = TRUE,   qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_4B_Reactome_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(Reactome1, showCategory = 20) )
  dev.off()
  
  png(filename=paste(MyFolder, "/", MyFileName, "_4B_Reactome_enrichment.png", sep=""),  width=800, height=500)
  print( dotplot(Reactome1, showCategory = 20) )
  dev.off()
  
  svg(filename=paste(MyFolder, "/", MyFileName, "_4B_Reactome_enrichment.svg", sep=""),  width=8, height=5)
  print( dotplot(Reactome1, showCategory = 20) )
  dev.off()
  
  #####################
  Reactome2     <- as.data.frame( Reactome1  )
  
  if(nrow(Reactome2) > 1) {
    Reactome2$GeneRatio <- sapply(Reactome1$GeneRatio, function(x) eval(parse(text=x))) * 100
    Reactome2$pvalue    <- -log10(Reactome2$pvalue)
    Reactome2$p.adjust  <- -log10(Reactome2$p.adjust)
    Reactome2$qvalue    <- -log10(Reactome2$qvalue)
    
    Reactome3 <- Reactome2[1:20, -8]
    if( nrow(Reactome2) <20 ) {Reactome3 <- Reactome2[, -8]}
    
    Reactome3 <- transform(Reactome3, ID= rev( factor(ID, levels=unique(ID))) )
    myMidValue <- ( max(Reactome3$GeneRatio) + min(Reactome3$GeneRatio)  )/2
    FigureTemp1 <- ggplot( data = Reactome3, aes(x = p.adjust, y = ID, size=Count , colour=GeneRatio ) ) + 
      geom_point( ) +  scale_colour_gradient2( low = "blue", mid = "yellow2", high = "red" ,  midpoint = myMidValue   ) +
      scale_y_discrete( breaks=Reactome3$ID, labels=rev(Reactome3$ID) ) + 
      xlab("-log10(ajusted p-value)") +   ylab("Pathway terms (Reactome)") +   ggtitle("Pathway terms (Reactome)") + 
      theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste( MyFileName, "_4C_Reactome_enrichment", sep=""),  height1=6,  width1=5)
  }
  
  
}
######################################################################################################################################################




sink( file=paste(outPath2,  "1-number-of-selected-peaks.txt",  sep="/")  )
######################################################################################################################################################
my_folder_5 <- paste(outPath2,  "5-peaks-annotation-all",  sep="/") 
if( ! file.exists(my_folder_5)  ) { dir.create(my_folder_5,  recursive = TRUE)  }
MyGeneOntology_1(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5,  MyFileName="5_allNearestGenes")

my_folder_6 <- paste(outPath2,  "6-peaks-annotation-100kb",  sep="/") 
if( ! file.exists(my_folder_6)  ) { dir.create(my_folder_6,  recursive = TRUE)  }
MyGeneOntology_1(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6,  MyFileName="6_100kb-NearestGenes")

my_folder_7 <- paste(outPath2,  "7-peaks-annotation-50kb",  sep="/") 
if( ! file.exists(my_folder_7)  ) { dir.create(my_folder_7,  recursive = TRUE)  }
MyGeneOntology_1(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7,  MyFileName="7_50kb-NearestGenes")

my_folder_8 <- paste(outPath2,  "8-peaks-annotation-10kb",  sep="/") 
if( ! file.exists(my_folder_8)  ) { dir.create(my_folder_8,  recursive = TRUE)  }
MyGeneOntology_1(MyMaxDis=10000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8,  MyFileName="8_10kb-NearestGenes")

my_folder_9 <- paste(outPath2,  "9-peaks-annotation-5kb",  sep="/") 
if( ! file.exists(my_folder_9)  ) { dir.create(my_folder_9,  recursive = TRUE)  }
MyGeneOntology_1(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9,  MyFileName="9_5kb-NearestGenes")

my_folder_10 <- paste(outPath2,  "10-peaks-annotation-1kb",  sep="/") 
if( ! file.exists(my_folder_10)  ) { dir.create(my_folder_10,  recursive = TRUE)  }
MyGeneOntology_1(MyMaxDis=1000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10,  MyFileName="10_1kb-NearestGenes")

######################################################################################################################################################
sink()





################################################################################
library(stringi)

sink( file=paste(outPath2,  "2-number-of-genes-forFurtherAnalysis.txt",  sep="/")  )

all_nearest_genes  <- as.data.frame(peakAnno)$geneId
all_nearest_genes2 <- do.call(paste, c(as.list(all_nearest_genes), sep=","))  
length(all_nearest_genes)
length(all_nearest_genes2)

nearest_genes_100kb  <- all_nearest_genes[ abs( as.data.frame(peakAnno)$distanceToTSS) < 100000 ]
nearest_genes_100kb_2 <- do.call(paste, c(as.list(nearest_genes_100kb), sep=","))  
length(nearest_genes_100kb)
length(nearest_genes_100kb_2)

nearest_genes_10kb  <- all_nearest_genes[ abs( as.data.frame(peakAnno)$distanceToTSS) < 10000 ]
nearest_genes_10kb_2 <- do.call(paste, c(as.list(nearest_genes_10kb), sep=","))  
length(nearest_genes_10kb)
length(nearest_genes_10kb_2)

sink()


shell_cmd1 <- paste("Rscript   clusterProfiler.oneList.hg38.R  ",   all_nearest_genes2 ,   "    ",  outPath2,   "/11_all_nearest_genes_moreAnalysis" ,   sep="") 
shell_cmd2 <- paste("Rscript   clusterProfiler.oneList.hg38.R  ",   nearest_genes_100kb_2, "    ",  outPath2,   "/12_100kb_nearest_genes_moreAnalysis",  sep="") 
shell_cmd3 <- paste("Rscript   clusterProfiler.oneList.hg38.R  ",   nearest_genes_10kb_2,  "    ",  outPath2,   "/13_10kb_nearest_genes_moreAnalysis",   sep="") 
system( shell_cmd1 )
system( shell_cmd2 )
system( shell_cmd3 )


 









