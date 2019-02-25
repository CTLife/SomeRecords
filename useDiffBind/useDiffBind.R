#! /usr/bin/env Rscript

##############################################################################################################################################################################################
suppressPackageStartupMessages( library(optparse) )  ## To run the script in command lines.

getParameters_f <- function() {
	option_list_Local <- list(   ## Options list with associated default value.  
  		optparse::make_option(opt_str=c("-I", "--inFile"),
			default="samplesInformation.example.txt",
			type="character",   dest="inFile",
			help="Path to the design/target file with samples information. [default: %default]."),

  		optparse::make_option(opt_str=c("-OD", "--outDir"),
			default="outDir_results",
			type="character",   dest="outDir",
                        help="Path to the directory containing all the analysis results. [default: %default]." )                                            
	)
	
	## Now parse the command line to check which option is given and get associated values.
	parser_Local <- optparse::OptionParser(usage="usage: %prog [options]",
			option_list=option_list_Local, 
			description="For DiffBind.",                             
			epilogue="For comments, bug reports etc..., please contact Yong Peng <yongp@outlook.com>."
	)
	opt_Local <- optparse::parse_args(parser_Local, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options
	return(opt_Local)
}
##############################################################################################################################################################################################





##############################################################################################################################################################################################
opt_g = getParameters_f()  

inFile_g  <- opt_g$inFile
outDir_g  <- opt_g$outDir


rm(getParameters_f)  
rm(opt_g)  

options(digits=10)  
suppressPackageStartupMessages( library(DiffBind) )
#  inFile_g = "H3K4me1_ICSIfresh-NC.csv"
#  outDir_g = "30-diffbind/All-H3K4me1/ICSIfresh-NC"
if( ! file.exists(outDir_g) ) { dir.create(outDir_g, recursive = TRUE) }  
##############################################################################################################################################################################################


samples_info <- read.csv( inFile_g )
write.table(x=samples_info, file = paste(outDir_g, "1-samples_info.txt", sep="/"),  quote = FALSE, 
            sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE )


DBA_object <- DiffBind::dba(sampleSheet=inFile_g)
sink( file = paste(outDir_g, "2a-DBA_object.txt", sep="/") )
   print(DBA_object)
sink()
pdf( file = paste(outDir_g, "2b-cross-correlations_occupancy.pdf", sep="/") )
    DiffBind::dba.plotHeatmap(DBA_object)  ## cross-correlations,  Correlation heatmap, using occupancy (peak caller score) data
    DiffBind::dba.plotPCA(DBA_object, attributes=c(DBA_CONDITION, DBA_TREATMENT) )  
dev.off()



binding_affinity_matrix <- DiffBind::dba.count(DBA_object)
sink( file = paste(outDir_g, "3a-binding_affinity_matrix.txt", sep="/") )
   print(binding_affinity_matrix)
sink()
pdf( file = paste(outDir_g, "3b-binding_affinity_matrix.pdf", sep="/") )
    DiffBind::dba.plotHeatmap(binding_affinity_matrix)  ## Correlation heatmap, using affinity (read count) data
    DiffBind::dba.plotPCA(binding_affinity_matrix, attributes=c(DBA_CONDITION, DBA_TREATMENT) )  
dev.off()


binding_affinity_contrast <- DiffBind::dba.contrast(binding_affinity_matrix, categories=DBA_CONDITION, minMembers=2)
sink( file = paste(outDir_g, "4-binding_affinity_contrast.txt", sep="/") )
print(binding_affinity_contrast)
sink()


binding_affinity_diff <- DiffBind::dba.analyze(binding_affinity_contrast, method=DBA_DESEQ2 )
sink( file = paste(outDir_g, "5a-binding_affinity_diff.txt", sep="/") )
print(binding_affinity_diff)
sink()
pdf( file = paste(outDir_g, "5b-binding_affinity_diff.pdf", sep="/") )
DiffBind::dba.plotHeatmap(binding_affinity_diff)   
DiffBind::dba.plotPCA(binding_affinity_diff, attributes=c(DBA_CONDITION, DBA_TREATMENT) )  
DiffBind::dba.plotHeatmap(binding_affinity_diff, contrast=1)  
DiffBind::dba.plotPCA(binding_affinity_diff, attributes=c(DBA_CONDITION, DBA_TREATMENT) , contrast=1)    
dev.off()





binding_affinity_diff.DB <- DiffBind::dba.report(binding_affinity_diff, th=1)

write.table(x=binding_affinity_diff.DB, file = paste(outDir_g, "6a-binding_affinity_diff.DB.all.txt", sep="/"),  quote = FALSE, 
            sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE )

binding_affinity_diff.DB2 <- DiffBind::dba.report(binding_affinity_diff, th=0.05)

write.table(x=binding_affinity_diff.DB2, file = paste(outDir_g, "6b-binding_affinity_diff.DB2.selected.txt", sep="/"),  quote = FALSE, 
            sep = "\t",  eol = "\n", na = "NA", dec = ".", row.names = FALSE,  col.names = TRUE )




pdf( file = paste(outDir_g, "7a-diff-figures.pdf", sep="/") )
DiffBind::dba.plotPCA(binding_affinity_diff,  label=DBA_ID, attributes=c(DBA_CONDITION, DBA_TREATMENT)  )   
DiffBind::dba.plotPCA(binding_affinity_diff,  label=DBA_ID, attributes=c(DBA_CONDITION, DBA_TREATMENT) , contrast=1)   
DiffBind::dba.plotMA(binding_affinity_diff)
DiffBind::dba.plotMA(binding_affinity_diff, bXY=TRUE)   
DiffBind::dba.plotVolcano(binding_affinity_diff)
DiffBind::dba.plotVolcano(binding_affinity_diff, th=0.05, fold=0.2630344 )  ## 2^fold
dev.off()

pdf( file = paste(outDir_g, "7b-diff-figures.pdf", sep="/") )
DiffBind::dba.plotBox(binding_affinity_diff)
DiffBind::dba.plotBox(binding_affinity_diff, attribute=DBA_TREATMENT)
DiffBind::dba.plotHeatmap(binding_affinity_diff, contrast=1, correlations=FALSE)
DiffBind::dba.plotHeatmap(binding_affinity_diff, contrast=1 )
dev.off()


png( file = paste(outDir_g, "7c-diff-figures.png", sep="/") )
DiffBind::dba.plotVolcano(binding_affinity_diff)
dev.off()

png( file = paste(outDir_g, "7d-diff-figures-1.2fold.png", sep="/") )
DiffBind::dba.plotVolcano(binding_affinity_diff, th=0.05, fold=0.2630344 )  ## 2^fold
dev.off()


olap.rate <- DiffBind::dba.overlap(binding_affinity_diff,mode=DBA_OLAP_RATE)
sink( file = paste(outDir_g, "8a-olap.rate.txt", sep="/") )
    print(olap.rate)
sink()
pdf( file = paste(outDir_g, "8b-olap.rate.pdf", sep="/") )
	plot(olap.rate, type='b',ylab='# peaks', xlab='Overlap at least this many peaksets')   
dev.off()



















