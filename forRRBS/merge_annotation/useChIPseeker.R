######################################################################################################################################################
##  The format of input file:
#1. The table header must be kept. The suffix of input files must be ".txt".
#2. The columns 1, 2 and 3 are required: chromosome, start, end. 
#3. The column 4: strand, (+, -, *)
#4. The column 5: name of DNA regions
#5. The column 6 "strength": peak strength, (1 or any real number)
#6. Such as "chr	start	end	strand	name	strength	pvalue	qvalue	meth.diff"
## Example: Rscript  useChIPseeker.R    inputDir    outDir   organism
######################################################################################################################################################





######################################################################################################################################################
args_g <- commandArgs(TRUE)
print("args: ")
print(args_g[1])         
print(args_g[2])
print(args_g[3])
print("#############")

inputDir   = args_g[1];     ## files of the same group in this folder (input dir)
outPath    = args_g[2];     ## output file path
RefGenome  = args_g[3];   
# inputDir = "5_DMR_Annotation/3-Formatted/forChIPseeker"
# outPath  = "5_DMR_Annotation/5-Annotation" 
# RefGenome  = "hg38"
if( ! file.exists(outPath)  ) { dir.create(outPath,  recursive = TRUE)  }
######################################################################################################################################################





######################################################################################################################################################
continue_on_error <- function() {
  print(" NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set 'options(error=continue_on_error())'. ")
}
# This is the key option
options(error=continue_on_error) 

suppressPackageStartupMessages( library(GenomicFeatures) ) 
suppressPackageStartupMessages( library(GenomicRanges) ) 
suppressPackageStartupMessages( library(clusterProfiler) ) 
suppressPackageStartupMessages( library(ReactomePA) ) 
suppressPackageStartupMessages( library(ChIPseeker) ) 
suppressPackageStartupMessages( library(DOSE) ) 
suppressPackageStartupMessages( library(ggplot2) ) 
suppressPackageStartupMessages( library(topGO) ) 
suppressPackageStartupMessages( library(KEGG.db) ) 
my_txdb_g  <- ""
my_orgdb_g <- ""
my_organism_g <- ""
my_organism2_g <- ""
if(RefGenome == "hg38") {
    suppressPackageStartupMessages( library(TxDb.Hsapiens.UCSC.hg38.knownGene) ) 
    suppressPackageStartupMessages( library(org.Hs.eg.db) )  
    my_txdb_g  <- TxDb.Hsapiens.UCSC.hg38.knownGene
    my_orgdb_g <- "org.Hs.eg.db"
    my_organism_g <- "human"
    my_organism2_g <- "hsa"
}
if(RefGenome == "mm10") {
    suppressPackageStartupMessages( library(TxDb.Mmusculus.UCSC.mm10.knownGene) ) 
    suppressPackageStartupMessages( library(org.Mm.eg.db) )  
    my_txdb_g  <- TxDb.Mmusculus.UCSC.mm10.knownGene
    my_orgdb_g <- "org.Mm.eg.db"
    my_organism_g <- "mouse"
    my_organism2_g <- "mmu"
}

sink( paste(outPath, "Parameters.txt", sep="/")   )
print( "my_txdb_g:"  )
print(  my_txdb_g  )
cat("\n\n\n\n\n") 
print( "my_orgdb_g: ")  
print(  my_orgdb_g ) 
cat("\n\n\n\n\n") 
print( "my_organism_g:" )
print( my_organism_g )
cat("\n\n\n\n\n") 
print( "my_organism2_g:" )
print( my_organism2_g )
cat("\n\n\n\n\n") 
sink()   

## read the input files.
peakFiles <- list.files(path = inputDir, pattern = ".txt",  full.names = TRUE )
peakFiles
length(peakFiles)
######################################################################################################################################################

 



######################################################################################################################################################
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



#############  Start the annotation for genomic regions.
MyPeaksAnno_OneGroup_1_g <- function(myPeak_t1, myPath_t1) { 
  myTempFunction <- function() {
  if( ! file.exists(myPath_t1)  ) { dir.create(myPath_t1,  recursive = TRUE)  }
  pdf(file=paste(myPath_t1, "1_ChIP_Peaks_over_Chromosomes.pdf", sep="/"),  width=20, height=20)
  print( covplot( myPeak_t1, weightCol="strength" , lower = 0.01  ) )
  dev.off()
  pdf(file=paste(myPath_t1, "2_ChIP_Peaks_over_Chromosomes.only-chromXY.pdf", sep="/"),  width=20, height=10)
  print( covplot( myPeak_t1, weightCol="strength" , chrs=c("chrX", "chrY") , lower = 0.01    ) )
  dev.off()
  pdf(file=paste(myPath_t1, "3_ChIP_Peaks_over_Chromosomes.chromLength.pdf", sep="/"),  width=20, height=20)
  print( covplot( myPeak_t1, weightCol="strength" , lower = -1  ) )
  dev.off()
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"MyPeaksAnno_OneGroup_1_g_111"}
  )
}


MyPeaksSignal_OneGroup_1_g <- function(myPeak_t1, myPath_t1, halfSize) { 
  myTempFunction <- function() {
  if( ! file.exists(myPath_t1)  ) { dir.create(myPath_t1,  recursive = TRUE)  }
  promoter2  <- getPromoters(TxDb=my_txdb_g, upstream=halfSize, downstream=halfSize)
  tagMatrix1 <- getTagMatrix(myPeak_t1, windows=promoter2)
  write.table(x=tagMatrix1, file =paste(myPath_t1, "1_peakMatrix_TSS.txt", sep="/"), 
              append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", 
              row.names = TRUE,  col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")

  pdf(file = paste(myPath_t1, "2A_Heatmap_of_ChIP_peaks_binding_to_TSS_regions.pdf", sep="/"),  width=10, height=20)
  print( tagHeatmap(tagMatrix1, xlim=c(-halfSize, halfSize), color="red") )
  dev.off()

    pdf(file=paste(myPath_t1, "2B_Average_Profile_of_ChIP_peaks_binding_to_TSS_region.pdf", sep="/"),  width=12, height=10)
  print( plotAvgProf(tagMatrix1, xlim=c(-halfSize, halfSize),  xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency" )  )
  dev.off()
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"MyPeaksSignal_OneGroup_1_g_222"}
  )
}


MyPeaksDistribution_OneGroup_1_g <- function(myPeak_t1, myPath_t1, halfSize=3000) { 
  myTempFunction <- function() {
    
  if( ! file.exists(myPath_t1)  ) { dir.create(myPath_t1,  recursive = TRUE)  }
  
  ## ------------------------------------------------------------------------
  peakAnno <- annotatePeak( myPeak_t1, tssRegion=c(-halfSize, halfSize),  TxDb=my_txdb_g, annoDb=my_orgdb_g  )  
  sink(  file = paste(myPath_t1, "1A_peakAnno_summary.txt", sep="/") )
  print( peakAnno )  
  sink()
  write.table(x=as.data.frame(peakAnno), file = paste(myPath_t1, "1B_peakAnno.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  ## ----fig.cap="Genomic Annotation by pieplot", fig.align="center", fig.height=6, fig.width=8----
  pdf(file=paste(myPath_t1, "2A_Genomic_Annotation_by_pieplot.pdf", sep="/"),  width=7, height=5)
  print( plotAnnoPie(peakAnno) ) 
  dev.off()
  
  svg(filename=paste(myPath_t1, "2B_Genomic_Annotation_by_pieplot.svg", sep="/"),  width=7, height=5)
  print( plotAnnoPie(peakAnno) ) 
  dev.off()
  
  ## ----fig.cap="Genomic Annotation by barplot", fig.align="center", fig.height=4, fig.width=10----
  pdf(file=paste(myPath_t1, "3A_Genomic_Annotation_by_barplot.pdf", sep="/"),  width=7, height=3)
  print( plotAnnoBar(peakAnno) ) 
  dev.off()
  
  svg(filename=paste(myPath_t1, "3B_Genomic_Annotation_by_barplot.svg", sep="/"),  width=7, height=3)
  print( plotAnnoBar(peakAnno) ) 
  dev.off()
  
  ## ----fig.cap="Genomic Annotation by vennpie", fig.align="center", fig.height=8, fig.width=11----
  pdf(file=paste(myPath_t1, "4A_Genomic_Annotation_by_vennpie.pdf", sep="/"),  width=7, height=5)
  print( vennpie(peakAnno) ) 
  dev.off()
  
  svg(filename=paste(myPath_t1, "4B_Genomic_Annotation_by_vennpie.svg", sep="/"),  width=7, height=5)
  print( vennpie(peakAnno) ) 
  dev.off()
  
  ## ----eval=F, fig.cap="Genomic Annotation by upsetplot", fig.align="center", fig.height=8, fig.width=12----
  pdf(file=paste(myPath_t1, "5A_Genomic_Annotation_by_upsetplot.pdf", sep="/"),  width=8, height=5)
  print( upsetplot(peakAnno) ) 
  dev.off()
  
  svg(filename=paste(myPath_t1, "5B_Genomic_Annotation_by_upsetplot.svg", sep="/"),  width=8, height=5)
  print( upsetplot(peakAnno) ) 
  dev.off()
  
  ## ----eval=F, fig.cap="Genomic Annotation by upsetplot", fig.align="center", fig.height=8, fig.width=12----
  pdf(file=paste(myPath_t1, "6A_Genomic_Annotation_by_upsetplot_vennpie.pdf", sep="/"),  width=8, height=5)
  print( upsetplot(peakAnno, vennpie=TRUE) ) 
  dev.off()
  
  svg(filename=paste(myPath_t1, "6B_Genomic_Annotation_by_upsetplot_vennpie.svg", sep="/"),  width=8, height=5)
  print( upsetplot(peakAnno, vennpie=TRUE) ) 
  dev.off()
  
  ## ----fig.cap="Distribution of Binding Sites", fig.align="center", fig.height=2, fig.width=6----
  pdf(file=paste(myPath_t1, "7A_Distribution_of_peaks_to_TSS.pdf", sep="/"),  width=8, height=2)
  print( plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS") ) 
  dev.off()
  
  svg(filename=paste(myPath_t1, "7B_Distribution_of_peaks_to_TSS.svg", sep="/"),  width=800, height=200)
  print( plotDistToTSS(peakAnno,  title="Distribution of peaks relative to TSS") ) 
  dev.off()
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"MyPeaksDistribution_OneGroup_1_g_333"}
  )
}


############# over-representation test (Enrichment Analysis)
MyGeneOntology_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction <- function() {
    
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS)<MyMaxDis, ]
  write.table(x=as.data.frame(peakAnno2), file = paste(MyFolder, "genomic_regions_annotation.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",       
              eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  print(nrow(peakAnno2))  ## number of peaks are selecred.
  
  bp1 <- enrichGO(gene=as.data.frame(peakAnno2)$geneId, OrgDb=my_orgdb_g, keyType = "ENTREZID", ont = "BP",
           pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2,
           minGSSize = 5, maxGSSize = 5000, readable = TRUE, pool = FALSE)
  emapplot(bp1)
  write.table(x= bp1 , file = paste(MyFolder, "/", MyFileName, "_1A_BP_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE,   qmethod = c("escape", "double"), fileEncoding = "")
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_1B_BP_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(bp1, showCategory = 20) )
  print( emapplot(bp1, showCategory = 20) )
  print( cnetplot(bp1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(bp1, showCategory = 20) )
  print( dotplot(bp1, showCategory = 10) )
  print( emapplot(bp1, showCategory = 10) )
  print( cnetplot(bp1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(bp1, showCategory = 10) )
  print( dotplot(bp1, showCategory = 5) )
  print( emapplot(bp1, showCategory = 5) )
  print( cnetplot(bp1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(bp1, showCategory = 5) )
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
  MF1 <- enrichGO(gene=as.data.frame(peakAnno2)$geneId, OrgDb=my_orgdb_g, keyType = "ENTREZID", ont = "MF",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2,
                  minGSSize = 5, maxGSSize = 5000, readable = TRUE, pool = FALSE)
  
  write.table(x= MF1 , file = paste(MyFolder, "/", MyFileName, "_2A_MF_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE,  col.names = TRUE,   qmethod = c("escape", "double"),  fileEncoding = "")
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_2B_MF_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(MF1, showCategory = 20) )
  print( emapplot(MF1, showCategory = 20) )
  print( cnetplot(MF1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(MF1, showCategory = 20) )
  print( dotplot(MF1, showCategory = 10) )
  print( emapplot(MF1, showCategory = 10) )
  print( cnetplot(MF1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(MF1, showCategory = 10) )
  print( dotplot(MF1, showCategory = 5) )
  print( emapplot(MF1, showCategory = 5) )
  print( cnetplot(MF1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(MF1, showCategory = 5) )
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
  CC1 <- enrichGO(gene=as.data.frame(peakAnno2)$geneId, OrgDb=my_orgdb_g, keyType = "ENTREZID", ont = "CC",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=NULL, qvalueCutoff = 0.2,
                  minGSSize = 5, maxGSSize = 5000, readable = TRUE, pool = FALSE)
  
  write.table(x= CC1 , file = paste(MyFolder, "/", MyFileName, "_3A_CC_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE,  col.names = TRUE,   qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_3B_CC_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(CC1, showCategory = 20) )
  print( emapplot(CC1, showCategory = 20) )
  print( cnetplot(CC1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(CC1, showCategory = 20) )
  print( dotplot(CC1, showCategory = 10) )
  print( emapplot(CC1, showCategory = 10) )
  print( cnetplot(CC1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(CC1, showCategory = 10) )
  print( dotplot(CC1, showCategory = 5) )
  print( emapplot(CC1, showCategory = 5) )
  print( cnetplot(CC1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
  print( goplot(CC1, showCategory = 5) ) 
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
  Reactome1 <- enrichPathway(gene=as.data.frame(peakAnno2)$geneId, organism = my_organism_g, pvalueCutoff = 0.05,
                pAdjustMethod = "BH", qvalueCutoff = 0.2, universe=NULL, minGSSize = 5, maxGSSize = 5000, readable = TRUE)
  
  write.table(x= Reactome1 , file = paste(MyFolder, "/", MyFileName, "_4A_Reactome_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE,  col.names = TRUE,   qmethod = c("escape", "double"),  fileEncoding = "")
  
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_4B_Reactome_enrichment.pdf", sep=""),  width=8, height=5)
  print( dotplot(Reactome1, showCategory = 20) )
  print( emapplot(Reactome1, showCategory = 20) )
  print( cnetplot(Reactome1, showCategory = 20, categorySize="qvalue", foldChange="GeneRatio") )
  print( dotplot(Reactome1, showCategory = 10) )
  print( emapplot(Reactome1, showCategory = 10) )
  print( cnetplot(Reactome1, showCategory = 10, categorySize="qvalue", foldChange="GeneRatio") )
  print( dotplot(Reactome1, showCategory = 5) )
  print( emapplot(Reactome1, showCategory = 5) )
  print( cnetplot(Reactome1, showCategory = 5, categorySize="qvalue", foldChange="GeneRatio") )
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
  tryCatch(
    myTempFunction(),
    error = function(err){"MyGeneOntology_OneGroup_1_g_444"}
  )
  
}


MyKEGG_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction <- function() {
    
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS)<MyMaxDis, ]
  print(nrow(peakAnno2))  ## number of peaks are selecred.
  
  enrich_KEGG <- enrichKEGG(gene=as.data.frame(peakAnno2)$geneId,  organism = my_organism2_g,  keyType = "kegg", 
                            pvalueCutoff = 0.05,  pAdjustMethod = "BH", universe=NULL, minGSSize = 5, maxGSSize = 5000,
                            qvalueCutoff = 0.2,   use_internal_data = TRUE )
  
  write.table(x= enrich_KEGG , file = paste(MyFolder, "/", MyFileName, "_1_KEGG_enrichment.txt", sep=""), append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n",  na = "NA",   dec = ".",   row.names = FALSE, col.names = TRUE,   qmethod = c("escape", "double"), fileEncoding = "")
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_2_KEGG_enrichment.pdf", sep=""),  width=12, height=10)
    print( barplot(enrich_KEGG, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "")  )
    print( barplot(enrich_KEGG, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "")  )
    print( dotplot(enrich_KEGG, x="geneRatio", color="p.adjust", showCategory=20, split=NULL, font.size=12, title="")  )
    print( dotplot(enrich_KEGG, x="count",     color="p.adjust", showCategory=20, split=NULL, font.size=12, title="")  )
    print( emapplot(enrich_KEGG, showCategory = 10, color = "p.adjust")  )
    print( emapplot(enrich_KEGG, showCategory = 20, color = "p.adjust")  )
    print( cnetplot(enrich_KEGG, showCategory = 2,   categorySize="pvalue", fixed = TRUE )  )
    print( cnetplot(enrich_KEGG, showCategory = 5,   categorySize="pvalue", fixed = TRUE )  )
    print( cnetplot(enrich_KEGG, showCategory = 10,  categorySize="pvalue", fixed = TRUE )  )
  dev.off()
  
  

  ##################################
  bp2     <- as.data.frame( enrich_KEGG ) 

  if(nrow(bp2) > 1) {
    bp2$GeneRatio <- sapply(enrich_KEGG$GeneRatio, function(x) eval(parse(text=x))) * 100
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
      xlab("-log10(ajusted p-value)") +   ylab("KEGG terms") +   ggtitle("KEGG terms") +  theme_bw()
    MySaveGgplot2_1( path1=MyFolder, fileName1=paste(MyFileName, "_1_KEGG_enrichment", sep=""),  height1=6,  width1=5)
  }
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"MyKEGG_OneGroup_1_g_555"}
  )
}


Disease_OneGroup_1_g <- function(MyMaxDis=900000000,  MyPeaksAnno, MyFolder,  MyFileName ) {
  myTempFunction <- function() {
    
  peakAnno2 <- as.data.frame(MyPeaksAnno)
  peakAnno2 <- peakAnno2[ abs(peakAnno2$distanceToTSS)<MyMaxDis, ]
  print(nrow(peakAnno2))  ## number of peaks are selecred.
  inputGenes2 = as.data.frame(peakAnno2)$geneId
  
  enrich_DO <- enrichDO(gene=inputGenes2, ont = "DO", pvalueCutoff = 0.05, pAdjustMethod = "BH",
           universe=NULL, minGSSize = 5, maxGSSize = 5000, qvalueCutoff = 0.2, readable = TRUE)
  
  write.table(x= enrich_DO , file = paste(MyFolder, "/", MyFileName, "_1_DO_enrichment.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
  
  pdf(file=paste(MyFolder, "/", MyFileName, "_2_DO_enrichment.pdf", sep=""),  width=12, height=10)
      print( barplot(enrich_DO, drop=TRUE, showCategory=20, x = "Count",     color = "p.adjust", font.size = 12, title = "") )
      print( barplot(enrich_DO, drop=TRUE, showCategory=20, x = "GeneRatio", color = "p.adjust", font.size = 12, title = "") )
      print( dotplot(enrich_DO, x="geneRatio", color="p.adjust", showCategory=20, split=NULL, font.size=12, title="") )
      print( dotplot(enrich_DO, x="count",     color="p.adjust", showCategory=20, split=NULL, font.size=12, title="") )
      print( emapplot(enrich_DO, showCategory = 10, color = "p.adjust")  )
      print( emapplot(enrich_DO, showCategory = 20, color = "p.adjust")  )
      print( cnetplot(enrich_DO, showCategory = 2,  categorySize="pvalue", fixed = TRUE ) )
      print( cnetplot(enrich_DO, showCategory = 5,  categorySize="pvalue", fixed = TRUE ) )
  dev.off()
  
  
  
  enrich_NCG <- enrichNCG(gene=inputGenes2,  pvalueCutoff=0.05, pAdjustMethod="BH", minGSSize=5, maxGSSize=500, qvalueCutoff=0.2,  readable=TRUE)
  write.table(x= enrich_NCG , file = paste(MyFolder, "/", MyFileName, "_3_NCG_enrichment.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
  pdf(file=paste(MyFolder, "/", MyFileName, "_4_NCG_enrichment.pdf", sep=""),  width=12, height=10)
    print( barplot(enrich_NCG, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "") )
    print( barplot(enrich_NCG, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "") )
    print( dotplot(enrich_NCG, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="") )
    print( dotplot(enrich_NCG, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="") )
    print( enrichMap(enrich_NCG, n = 10, fixed = TRUE, vertex.label.font = 1) )
    print( enrichMap(enrich_NCG, n = 20, fixed = TRUE, vertex.label.font = 1) )
    print( cnetplot(enrich_NCG, showCategory = 2,  categorySize="pvalue", fixed = TRUE ) )
    print( cnetplot(enrich_NCG, showCategory = 5,  categorySize="pvalue", fixed = TRUE ) )
  dev.off()
  
  
  
  
  enrich_DGN <- enrichDGN(gene=inputGenes2,  pvalueCutoff=0.05, pAdjustMethod="BH", minGSSize=5, maxGSSize=500, qvalueCutoff=0.2,  readable=TRUE)
  write.table(x= enrich_DGN , file = paste(MyFolder, "/", MyFileName, "_5_DGN_enrichment.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
  pdf(file=paste(MyFolder, "/", MyFileName, "_6_DGN_enrichment.pdf", sep=""),  width=12, height=10)  
      print( barplot(enrich_DGN, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "") )
      print( barplot(enrich_DGN, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "") )
      print( dotplot(enrich_DGN, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="") )
      print( dotplot(enrich_DGN, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="") )
      print( enrichMap(enrich_DGN, n = 10, fixed = TRUE, vertex.label.font = 1) )
      print( enrichMap(enrich_DGN, n = 20, fixed = TRUE, vertex.label.font = 1) )
      print( cnetplot(enrich_DGN, showCategory = 2,  categorySize="pvalue", fixed = TRUE ) )
      print( cnetplot(enrich_DGN, showCategory = 5,  categorySize="pvalue", fixed = TRUE ) )
  dev.off()
  }
  tryCatch(
    myTempFunction(),
    error = function(err){"Disease_OneGroup_1_g_666"}
  )
}




######################################################################################################################################################




MyAnnotation_OneGroup_1_g <- function(myFile_t1, myPath_t1) {
  ## myFile_t1 = peakFiles[1] 
  ## myPath_t1 = "3-annotation/group_1"
  if( ! file.exists( myPath_t1 )  ) { dir.create(myPath_t1,  recursive = TRUE)  }
  sink(file = paste(myPath_t1, "1_fileName.txt", sep="/") )
  print( myFile_t1 )
  sink()
  myPeak_t1 = readPeakFile( myFile_t1 )
  
  MyPeaksAnno_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "1-peaks-on-chromosomes",  sep="/") )
  MyPeaksSignal_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "2A-peaks-on-TSSs-10kb",  sep="/"), halfSize = 10000)
  #MyPeaksSignal_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "2B-peaks-on-TSSs-4kb",  sep="/"), halfSize = 4000)
  MyPeaksDistribution_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "3A-peaks-distribution-TSS5kb",  sep="/"), halfSize = 5000)
  MyPeaksDistribution_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "3B-peaks-distribution-TSS3kb",  sep="/"), halfSize = 3000)
  MyPeaksDistribution_OneGroup_1_g(myPeak_t1 = myPeak_t1, myPath_t1 = paste(myPath_t1,  "3C-peaks-distribution-TSS1kb",  sep="/"), halfSize = 1000)
  
  peakAnno <- annotatePeak( myPeak_t1, tssRegion=c(-1000, 1000),  TxDb=my_txdb_g, annoDb=my_orgdb_g  )  
  write.table(x=as.data.frame(peakAnno), file = paste(myPath_t1, "genomic_regions_annotation_TSS1kb.txt", sep="/"), append = FALSE, quote = FALSE, sep = "\t",       
              eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"),  fileEncoding = "")
  
  sink( file=paste(myPath_t1,  "2-number-of-selected-peaks.txt",  sep="/")  )
  
  my_folder_5 <- paste(myPath_t1,  "4A-GeneOntology-all",  sep="/") 
  if( ! file.exists(my_folder_5)  ) { dir.create(my_folder_5,  recursive = TRUE)  }
  MyGeneOntology_OneGroup_1_g(MyMaxDis=900000000, MyPeaksAnno=peakAnno, MyFolder=my_folder_5,  MyFileName="4A_allNearestGenes")
  
  my_folder_6 <- paste(myPath_t1,  "4B-GeneOntology-100kb",  sep="/") 
  if( ! file.exists(my_folder_6)  ) { dir.create(my_folder_6,  recursive = TRUE)  }
  MyGeneOntology_OneGroup_1_g(MyMaxDis=100000, MyPeaksAnno=peakAnno, MyFolder=my_folder_6,  MyFileName="4B_100kb-NearestGenes")
  
  my_folder_7 <- paste(myPath_t1,  "4C-GeneOntology-50kb",  sep="/") 
  if( ! file.exists(my_folder_7)  ) { dir.create(my_folder_7,  recursive = TRUE)  }
  MyGeneOntology_OneGroup_1_g(MyMaxDis=50000, MyPeaksAnno=peakAnno, MyFolder=my_folder_7,  MyFileName="4C_50kb-NearestGenes")
  
  my_folder_8 <- paste(myPath_t1,  "4D-GeneOntology-10kb",  sep="/") 
  if( ! file.exists(my_folder_8)  ) { dir.create(my_folder_8,  recursive = TRUE)  }
  MyGeneOntology_OneGroup_1_g(MyMaxDis=10000, MyPeaksAnno=peakAnno, MyFolder=my_folder_8,  MyFileName="4D_10kb-NearestGenes")
  
  my_folder_9 <- paste(myPath_t1,  "4E-GeneOntology-5kb",  sep="/") 
  if( ! file.exists(my_folder_9)  ) { dir.create(my_folder_9,  recursive = TRUE)  }
  MyGeneOntology_OneGroup_1_g(MyMaxDis=5000, MyPeaksAnno=peakAnno, MyFolder=my_folder_9,  MyFileName="4E_5kb-NearestGenes")
  
  my_folder_10 <- paste(myPath_t1,  "4F-GeneOntology-20kb",  sep="/") 
  if( ! file.exists(my_folder_10)  ) { dir.create(my_folder_10,  recursive = TRUE)  }
  MyGeneOntology_OneGroup_1_g(MyMaxDis=20000, MyPeaksAnno=peakAnno, MyFolder=my_folder_10,  MyFileName="4F_20kb-NearestGenes")
  
  sink()
  
  
  my_folder_11 <- paste(myPath_t1,  "5A-KEGG_all",  sep="/") 
  if( ! file.exists(my_folder_11)  ) { dir.create(my_folder_11,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=900000000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_11,  MyFileName="5A_allNearestGenes" )
    
  
  my_folder_12 <- paste(myPath_t1,  "5B-KEGG_100kb",  sep="/") 
  if( ! file.exists(my_folder_12)  ) { dir.create(my_folder_12,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=100000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_12,  MyFileName="5B_100kb-NearestGenes" )
  
  
  my_folder_13 <- paste(myPath_t1,  "5C-KEGG_50kb",  sep="/") 
  if( ! file.exists(my_folder_13)  ) { dir.create(my_folder_13,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=50000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_13,  MyFileName="5C_50kb-NearestGenes" )
  
  
  my_folder_14 <- paste(myPath_t1,  "5D-KEGG_10kb",  sep="/") 
  if( ! file.exists(my_folder_14)  ) { dir.create(my_folder_14,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=10000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_14,  MyFileName="5D_10kb-NearestGenes" )
  
  
  my_folder_15 <- paste(myPath_t1,  "5E-KEGG_5kb",  sep="/") 
  if( ! file.exists(my_folder_15)  ) { dir.create(my_folder_15,  recursive = TRUE)  }
  MyKEGG_OneGroup_1_g(MyMaxDis=10000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_15,  MyFileName="5E_5kb-NearestGenes" )
  
  
  
  my_folder_16 <- paste(myPath_t1,  "6A-Disease_all",  sep="/") 
  if( ! file.exists(my_folder_16)  ) { dir.create(my_folder_16,  recursive = TRUE)  }
  Disease_OneGroup_1_g(MyMaxDis=900000000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_16,  MyFileName="6A_allNearestGenes" )
  
  my_folder_16 <- paste(myPath_t1,  "6B-Disease_100kb",  sep="/") 
  if( ! file.exists(my_folder_16)  ) { dir.create(my_folder_16,  recursive = TRUE)  }
  Disease_OneGroup_1_g(MyMaxDis=100000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_16,  MyFileName="6B_100kb-NearestGenes" )
  
  my_folder_16 <- paste(myPath_t1,  "6C-Disease_50kb",  sep="/") 
  if( ! file.exists(my_folder_16)  ) { dir.create(my_folder_16,  recursive = TRUE)  }
  Disease_OneGroup_1_g(MyMaxDis=50000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_16,  MyFileName="6C_50kb-NearestGenes" )
  
  my_folder_16 <- paste(myPath_t1,  "6D-Disease_10kb",  sep="/") 
  if( ! file.exists(my_folder_16)  ) { dir.create(my_folder_16,  recursive = TRUE)  }
  Disease_OneGroup_1_g(MyMaxDis=10000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_16,  MyFileName="6D_10kb-NearestGenes" )
  
  my_folder_16 <- paste(myPath_t1,  "6E-Disease_5kb",  sep="/") 
  if( ! file.exists(my_folder_16)  ) { dir.create(my_folder_16,  recursive = TRUE)  }
  Disease_OneGroup_1_g(MyMaxDis=5000,  MyPeaksAnno=peakAnno, MyFolder=my_folder_16,  MyFileName="6E_5kb-NearestGenes" )
  
  
}




for(i in c(1:length(peakFiles)) ) {
  # i = 1
  MyAnnotation_OneGroup_1_g( myFile_t1 = peakFiles[i],  myPath_t1 = paste(outPath, "/File_", i, sep="") )
}




















