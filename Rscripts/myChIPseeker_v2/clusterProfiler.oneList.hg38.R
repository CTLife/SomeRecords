## example: Rscript   clusterProfiler.oneList.hg38.R     geneList_vector    outputPath  

args <- commandArgs(TRUE)
print("args: ")
print(args[1])         
print(args[2])      
print("#############")

inputGenes   = args[1];  # ENTREZID (GeneID)
outputPath   = args[2];  # ouput path of all figures and results 
inputGenes   = strsplit(inputGenes, ",")
inputGenes   = as.vector(inputGenes[[1]])

# m1 <- read.table(file="4_peakAnno.txt", header = TRUE, sep = "\t", quote = "", dec = "."  )
# dim(m1)
# m1[1:20,]
# inputGenes   = m1[,15]
# inputGenes[1:20]
# length(inputGenes)
# outputPath   = "4-overlap-pups-annotation/ART-E18-H3K4me1/2_ART-E18D-H3K4me1.bed"

print("##################################### 1 #################################################")
outNumGenes  = paste("Number of Genes:", length(inputGenes) ,   sep="    ")  
print(outNumGenes)
if( ! file.exists(outputPath) ) { dir.create(outputPath, recursive = TRUE)  }

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(pathview)
library(ggplot2)
library(topGO)
my_TxDb <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
my_orgdb <- "org.Hs.eg.db"

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


inputGenes2 <- as.character(inputGenes)
inputGenes[1:20]
inputGenes2[1:20]




continue_on_error <- function()
{
  print("NOTE: THERE WAS AN ERROR HERE. We are continuing because we have set
'options(error=continue_on_error())'")
}

# This is the key option
options(error=continue_on_error) 



############################################## 1. GO Analysis
outputPath2 <- paste(outputPath, "/1-GO-Analysis",   sep="")  
if( ! file.exists(outputPath2) ) { dir.create(outputPath2, recursive = TRUE)  }
 
## 1.1 GO classification
GO_class_BP_1 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "BP", level = 1, readable = TRUE)
write.table(x= GO_class_BP_1 , file = paste(outputPath2, "/1.1_1_GO_class_BP_level-1.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_BP_2 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "BP", level = 2, readable = TRUE)
write.table(x= GO_class_BP_2 , file = paste(outputPath2, "/1.1_1_GO_class_BP_level-2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_BP_3 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "BP", level = 3, readable = TRUE)
write.table(x= GO_class_BP_3 , file = paste(outputPath2, "/1.1_1_GO_class_BP_level-3.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_BP_4 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "BP", level = 4, readable = TRUE)
write.table(x= GO_class_BP_4 , file = paste(outputPath2, "/1.1_1_GO_class_BP_level-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_BP_5 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "BP", level = 5, readable = TRUE)
write.table(x= GO_class_BP_5 , file = paste(outputPath2, "/1.1_1_GO_class_BP_level-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )


######
GO_class_MF_1 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "MF", level = 1, readable = TRUE)
write.table(x= GO_class_MF_1 , file = paste(outputPath2, "/1.1_2_GO_class_MF_level-1.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_MF_2 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "MF", level = 2, readable = TRUE)
write.table(x= GO_class_MF_2 , file = paste(outputPath2, "/1.1_2_GO_class_MF_level-2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_MF_3 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "MF", level = 3, readable = TRUE)
write.table(x= GO_class_MF_3 , file = paste(outputPath2, "/1.1_2_GO_class_MF_level-3.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_MF_4 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "MF", level = 4, readable = TRUE)
write.table(x= GO_class_MF_4 , file = paste(outputPath2, "/1.1_2_GO_class_MF_level-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_MF_5 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "MF", level = 5, readable = TRUE)
write.table(x= GO_class_MF_5 , file = paste(outputPath2, "/1.1_2_GO_class_MF_level-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )


######
GO_class_CC_1 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "CC", level = 1, readable = TRUE)
write.table(x= GO_class_CC_1 , file = paste(outputPath2, "/1.1_3_GO_class_CC_level-1.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_CC_2 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "CC", level = 2, readable = TRUE)
write.table(x= GO_class_CC_2 , file = paste(outputPath2, "/1.1_3_GO_class_CC_level-2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_CC_3 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "CC", level = 3, readable = TRUE)
write.table(x= GO_class_CC_3 , file = paste(outputPath2, "/1.1_3_GO_class_CC_level-3.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_CC_4 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "CC", level = 4, readable = TRUE)
write.table(x= GO_class_CC_4 , file = paste(outputPath2, "/1.1_3_GO_class_CC_level-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )

GO_class_CC_5 <- groupGO(gene = inputGenes2, OrgDb = my_orgdb, keytype="ENTREZID", ont = "CC", level = 5, readable = TRUE)
write.table(x= GO_class_CC_5 , file = paste(outputPath2, "/1.1_3_GO_class_CC_level-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )




## 1.2 GO over-representation test
enrichGO_BP <- enrichGO(gene=inputGenes2,  OrgDb=my_orgdb, keytype="ENTREZID", ont="BP", pvalueCutoff = 0.01,  pAdjustMethod = "BH",  # universe, 
                qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 1000,  readable=TRUE, pool=FALSE)
write.table(x= enrichGO_BP , file = paste(outputPath2, "/1.2_1_enrichGO_BP.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1_enrichGO_BP.pdf", sep=""),  width=12, height=10)
    barplot(enrichGO_BP, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
    barplot(enrichGO_BP, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
    dotplot(enrichGO_BP, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
    dotplot(enrichGO_BP, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
    enrichMap(enrichGO_BP, n = 10, fixed = TRUE, vertex.label.font = 1)
    enrichMap(enrichGO_BP, n = 20, fixed = TRUE, vertex.label.font = 1)
    cnetplot(enrichGO_BP, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
    cnetplot(enrichGO_BP, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
    plotGOgraph(enrichGO_BP,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
    plotGOgraph(enrichGO_BP,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_BP2 <- dropGO(enrichGO_BP, level = c(1:4), term = NULL)  ## drop specific GO terms or level
write.table(x= enrichGO_BP2 , file = paste(outputPath2, "/1.2_1A_enrichGO_BP_dropLevel1-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1A_enrichGO_BP_dropLevel1-4.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_BP2, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_BP2, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_BP2, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_BP2, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_BP2, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_BP2, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_BP2, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_BP2, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_BP2,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_BP2,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_BP3 <- dropGO(enrichGO_BP, level = c(1:5), term = NULL)  ## drop specific GO terms or level
write.table(x= enrichGO_BP3 , file = paste(outputPath2, "/1.2_1B_enrichGO_BP_dropLevel1-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1B_enrichGO_BP_dropLevel1-5.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_BP3, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_BP3, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_BP3, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_BP3, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_BP3, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_BP3, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_BP3, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_BP3, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_BP3,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_BP3,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_BP4 <- gofilter(enrichGO_BP, level=4)   ## test GO at sepcific level
write.table(x= enrichGO_BP4 , file = paste(outputPath2, "/1.2_1C_enrichGO_BP_Level-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1C_enrichGO_BP_Level-4.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_BP4, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_BP4, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_BP4, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_BP4, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_BP4, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_BP4, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_BP4, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_BP4, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_BP4,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_BP4,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()

 
enrichGO_BP5 <- gofilter(enrichGO_BP, level=5)   ## test GO at sepcific level
write.table(x= enrichGO_BP5 , file = paste(outputPath2, "/1.2_1D_enrichGO_BP_Level-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1D_enrichGO_BP_Level-5.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_BP5, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_BP5, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_BP5, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_BP5, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_BP5, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_BP5, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_BP5, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_BP5, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_BP5,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_BP5,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_BP6 <- simplify(enrichGO_BP, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)  ##  reduce redundancy of enriched GO terms
write.table(x= enrichGO_BP6 , file = paste(outputPath2, "/1.2_1E_enrichGO_BP_simplify.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1E_enrichGO_BP_simplify.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_BP6, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_BP6, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_BP6, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_BP6, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_BP6, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_BP6, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_BP6, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_BP6, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_BP6,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_BP6,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_BP7 <- simplify(enrichGO_BP, cutoff = 0.8, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)  ##  reduce redundancy of enriched GO terms
write.table(x= enrichGO_BP7 , file = paste(outputPath2, "/1.2_1F_enrichGO_BP_simplify2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_1F_enrichGO_BP_simplify2.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_BP7, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_BP7, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_BP7, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_BP7, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_BP7, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_BP7, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_BP7, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_BP7, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_BP7,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_BP7,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()



#################################
enrichGO_MF <- enrichGO(gene=inputGenes2,  OrgDb=my_orgdb, keytype="ENTREZID", ont="MF", pvalueCutoff = 0.01,  pAdjustMethod = "BH",  # universe, 
                        qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 1000,  readable=TRUE, pool=FALSE)
write.table(x= enrichGO_MF , file = paste(outputPath2, "/1.2_2_enrichGO_MF.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_2_enrichGO_MF.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_MF2 <- dropGO(enrichGO_MF, level = c(1:4), term = NULL)  ## drop specific GO terms or level
write.table(x= enrichGO_MF2 , file = paste(outputPath2, "/1.2_2A_enrichGO_MF_dropLevel1-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
  pdf(file=paste(outputPath2, "/1.2_2A_enrichGO_MF_dropLevel1-4.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF2, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF2, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF2, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF2, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF2, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF2, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF2, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF2, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF2,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF2,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()

enrichGO_MF3 <- dropGO(enrichGO_MF, level = c(1:5), term = NULL)  ## drop specific GO terms or level
write.table(x= enrichGO_MF3 , file = paste(outputPath2, "/1.2_2B_enrichGO_MF_dropLevel1-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
  pdf(file=paste(outputPath2, "/1.2_2B_enrichGO_MF_dropLevel1-5.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF3, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF3, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF3, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF3, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF3, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF3, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF3, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF3, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF3,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF3,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()



enrichGO_MF4 <- gofilter(enrichGO_MF, level=4)   ## test GO at sepcific level
write.table(x= enrichGO_MF4 , file = paste(outputPath2, "/1.2_2C_enrichGO_MF_Level-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_2C_enrichGO_MF_Level-4.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF4, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF4, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF4, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF4, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF4, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF4, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF4, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF4, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF4,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF4,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_MF5 <- gofilter(enrichGO_MF, level=5)   ## test GO at sepcific level
write.table(x= enrichGO_MF5 , file = paste(outputPath2, "/1.2_2D_enrichGO_MF_Level-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_2D_enrichGO_MF_Level-5.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF5, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF5, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF5, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF5, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF5, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF5, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF5, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF5, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF5,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF5,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_MF6 <- simplify(enrichGO_MF, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)  ##  reduce redundancy of enriched GO terms
write.table(x= enrichGO_MF6 , file = paste(outputPath2, "/1.2_2E_enrichGO_MF_simplify.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_2E_enrichGO_MF_simplify.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF6, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF6, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF6, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF6, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF6, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF6, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF6, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF6, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF6,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF6,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()



enrichGO_MF7 <- simplify(enrichGO_MF, cutoff = 0.8, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)  ##  reduce redundancy of enriched GO terms
write.table(x= enrichGO_MF7 , file = paste(outputPath2, "/1.2_2F_enrichGO_MF_simplify2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_2F_enrichGO_MF_simplify2.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_MF7, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_MF7, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_MF7, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_MF7, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_MF7, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_MF7, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_MF7, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_MF7, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_MF7,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_MF7,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()



#################################
enrichGO_CC <- enrichGO(gene=inputGenes2,  OrgDb=my_orgdb, keytype="ENTREZID", ont="CC", pvalueCutoff = 0.01,  pAdjustMethod = "BH",  # universe, 
                        qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 1000,  readable=TRUE, pool=FALSE)
write.table(x= enrichGO_CC , file = paste(outputPath2, "/1.2_3_enrichGO_CC.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3_enrichGO_CC.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_CC2 <- dropGO(enrichGO_CC, level = c(1:4), term = NULL)  ## drop specific GO terms or level
write.table(x= enrichGO_CC2 , file = paste(outputPath2, "/1.2_3A_enrichGO_CC_dropLevel1-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3A_enrichGO_CC_dropLevel1-4.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC2, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC2, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC2, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC2, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC2, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC2, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC2, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC2, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC2,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC2,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()

enrichGO_CC3 <- dropGO(enrichGO_CC, level = c(1:5), term = NULL)  ## drop specific GO terms or level
write.table(x= enrichGO_CC3 , file = paste(outputPath2, "/1.2_3B_enrichGO_CC_dropLevel1-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3B_enrichGO_CC_dropLevel1-5.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC3, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC3, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC3, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC3, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC3, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC3, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC3, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC3, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC3,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC3,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()



enrichGO_CC4 <- gofilter(enrichGO_CC, level=4)   ## test GO at sepcific level
write.table(x= enrichGO_CC4 , file = paste(outputPath2, "/1.2_3C_enrichGO_CC_Level-4.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3C_enrichGO_CC_Level-4.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC4, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC4, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC4, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC4, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC4, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC4, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC4, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC4, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC4,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC4,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_CC5 <- gofilter(enrichGO_CC, level=5)   ## test GO at sepcific level
write.table(x= enrichGO_CC5 , file = paste(outputPath2, "/1.2_3D_enrichGO_CC_Level-5.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3D_enrichGO_CC_Level-5.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC5, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC5, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC5, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC5, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC5, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC5, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC5, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC5, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC5,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC5,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


enrichGO_CC6 <- simplify(enrichGO_CC, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)  ##  reduce redundancy of enriched GO terms
write.table(x= enrichGO_CC6 , file = paste(outputPath2, "/1.2_3E_enrichGO_CC_simplify.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3E_enrichGO_CC_simplify.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC6, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC6, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC6, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC6, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC6, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC6, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC6, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC6, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC6,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC6,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()



enrichGO_CC7 <- simplify(enrichGO_CC, cutoff = 0.8, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)  ##  reduce redundancy of enriched GO terms
write.table(x= enrichGO_CC7 , file = paste(outputPath2, "/1.2_3F_enrichGO_CC_simplify2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath2, "/1.2_3F_enrichGO_CC_simplify2.pdf", sep=""),  width=12, height=10)
  barplot(enrichGO_CC7, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrichGO_CC7, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrichGO_CC7, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrichGO_CC7, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrichGO_CC7, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrichGO_CC7, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrichGO_CC7, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrichGO_CC7, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
  plotGOgraph(enrichGO_CC7,  firstSigNodes = 10,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
  plotGOgraph(enrichGO_CC7,  firstSigNodes = 20,  useInfo = "", sigForAll = TRUE, useFullNames = TRUE)
dev.off()


## 1.3  GO Gene Set Enrichment Analysis (GSEA), gene expression level is required.
# GO_GSEA_BP <- gseGO(geneList = inputGenes2,  OrgDb = my_orgdb, ont = "BP",  keyType = "ENTREZID", nPerm = 1000,
#              minGSSize = 100, maxGSSize = 1000, pvalueCutoff = 0.01,  by = "fgsea")  
# write.table(x= GO_GSEA_BP , file = paste(outputPath2, "/1.3_1_GO_GSEA_BP.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    


## 1.4 GO Semantic Similarity Analysis









############################################## 2. KEGG analysis
outputPath3 <- paste(outputPath, "/2-KEGG-Analysis",   sep="")  
if( ! file.exists(outputPath3) ) { dir.create(outputPath3, recursive = TRUE)  }

## 2.1 KEGG over-representation test
enrich_KEGG <- enrichKEGG(gene=inputGenes2, organism="hsa", keyType="kegg",  pvalueCutoff=0.05)
write.table(x= enrich_KEGG , file = paste(outputPath3, "/2.1_1_enrich_KEGG.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath3, "/2.1_1_enrich_KEGG.pdf", sep=""),  width=12, height=10)
  barplot(enrich_KEGG, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrich_KEGG, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrich_KEGG, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrich_KEGG, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrich_KEGG, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrich_KEGG, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrich_KEGG, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrich_KEGG, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
dev.off()


# browseKEGG(enrich_KEGG, "hsa04110")
# hsa04110 <- pathview(gene.data  = enrich_KEGG,  pathway.id = "hsa04110", species    = "hsa" )


## 2.2 KEGG Module over-representation test
enrich_KEGG2 <- enrichMKEGG(gene = inputGenes2, organism="hsa", keyType = "kegg",  pvalueCutoff = 0.05)
write.table(x= enrich_KEGG2 , file = paste(outputPath3, "/2.2_1_enrich_Module_KEGG2.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath3, "/2.2_1_enrich_Module_KEGG2.pdf", sep=""),  width=12, height=10)
  barplot(enrich_KEGG2, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrich_KEGG2, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrich_KEGG2, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrich_KEGG2, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrich_KEGG2, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrich_KEGG2, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrich_KEGG2, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrich_KEGG2, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
dev.off()









############################################## 3. Disease analysis
outputPath4 <- paste(outputPath, "/3-Disease-Analysis",   sep="")  
if( ! file.exists(outputPath4) ) { dir.create(outputPath4, recursive = TRUE)  }

enrich_DO <- enrichDO(gene=inputGenes2, ont="DO", pvalueCutoff=0.2, pAdjustMethod="BH", minGSSize=5, maxGSSize=500, qvalueCutoff=0.1,  readable=TRUE)
write.table(x= enrich_DO , file = paste(outputPath4, "/3.1_1_enrich_DO.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath4, "/3.1_1_enrich_DO.pdf", sep=""),  width=12, height=10)
  barplot(enrich_DO, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrich_DO, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrich_DO, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrich_DO, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrich_DO, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrich_DO, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrich_DO, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrich_DO, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
dev.off()



enrich_NCG <- enrichNCG(gene=inputGenes2,  pvalueCutoff=0.05, pAdjustMethod="BH", minGSSize=5, maxGSSize=500, qvalueCutoff=0.2,  readable=TRUE)
write.table(x= enrich_NCG , file = paste(outputPath4, "/3.2_1_enrich_NCG.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath4, "/3.2_1_enrich_NCG.pdf", sep=""),  width=12, height=10)
  barplot(enrich_NCG, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
  barplot(enrich_NCG, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
  dotplot(enrich_NCG, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  dotplot(enrich_NCG, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
  enrichMap(enrich_NCG, n = 10, fixed = TRUE, vertex.label.font = 1)
  enrichMap(enrich_NCG, n = 20, fixed = TRUE, vertex.label.font = 1)
  cnetplot(enrich_NCG, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
  cnetplot(enrich_NCG, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
dev.off()




enrich_DGN <- enrichDGN(gene=inputGenes2,  pvalueCutoff=0.05, pAdjustMethod="BH", minGSSize=5, maxGSSize=500, qvalueCutoff=0.2,  readable=TRUE)
write.table(x= enrich_DGN , file = paste(outputPath4, "/3.3_1_enrich_DGN.txt", sep=""),  append = FALSE, quote = FALSE, sep = "\t",  row.names = FALSE, col.names = TRUE )                                    
pdf(file=paste(outputPath4, "/3.3_1_enrich_DGN.pdf", sep=""),  width=12, height=10)
barplot(enrich_DGN, drop=TRUE, showCategory=20, x = "Count",     colorBy = "p.adjust", font.size = 12, title = "")
barplot(enrich_DGN, drop=TRUE, showCategory=20, x = "GeneRatio", colorBy = "p.adjust", font.size = 12, title = "")
dotplot(enrich_DGN, x="geneRatio", colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
dotplot(enrich_DGN, x="count",     colorBy="p.adjust", showCategory=20, split=NULL, font.size=12, title="")
enrichMap(enrich_DGN, n = 10, fixed = TRUE, vertex.label.font = 1)
enrichMap(enrich_DGN, n = 20, fixed = TRUE, vertex.label.font = 1)
cnetplot(enrich_DGN, showCategory = 2,  categorySize="pvalue", fixed = TRUE )
cnetplot(enrich_DGN, showCategory = 5,  categorySize="pvalue", fixed = TRUE )
dev.off()



 


