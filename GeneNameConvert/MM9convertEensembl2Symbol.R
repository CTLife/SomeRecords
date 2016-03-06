#! /usr/bin/env Rscript

## Convert ensembl id in the table to gene symbol in mm9 ( mmusculus_gene_ensembl Mus musculus genes (NCBIM37) NCBIM37 ).
## args[1]: file name, only csv or txt
## args[2]: column number of ensembl_id
## args[3]: if "1", input file with header; other: no header.

## for other reference genomes using: 
##                                   ensembl67=useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL')
##                                   listDatasets(ensembl67)[grep("Mus", listDatasets(ensembl67)$description),]
##                                   listDatasets(ensembl67)
## run this script, such as:  Rscript  MM9convertEensembl2Symbol.R   DoulbleKOvsWT.down.txt  1  1


args <- commandArgs(TRUE)

print("args: ")
print(args[1])
print(args[2])
print(args[3])
print("#############")

ensembl.col <- as.numeric(args[2])
print( ensembl.col )

if (args[3] == "1"){
	header = TRUE
}else{
	header = FALSE
}

## read file by suffix of file
require(tools)
file.type <- file_ext(args[1])

if (file.type == "txt"){
	  file.name   <- sub('.txt$', '', basename(args[1]))
	  input.table <- read.table(args[1], header=header, sep="\t", comment.char="", quote="", stringsAsFactors=FALSE)
}else if(file.type == "csv"){
	  file.name   <- sub('.csv$', '', basename(args[1]))
	  input.table <- read.csv(args[1], header=header, comment.char="", stringsAsFactors=FALSE)
}else{
	stop("Only csv and txt are supported now.")
}




require(biomaRt)
mart <- useMart(host='may2012.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',  dataset="mmusculus_gene_ensembl")
print( mart )


#input.table[, ensembl.col]

sink( paste(args[1], ".runLog.txt", sep="") )
listAttributes(mart)
listFilters(mart)
sink()
filterName <- "ensembl_transcript_id"   ## for others, by using listFilters(mart)


results <- getBM(attributes = c("ensembl_gene_id",  "ensembl_transcript_id",  "description", 
                                "mgi_symbol", "refseq_mrna", "refseq_ncrna", "external_gene_name" ),  
                 filters = c(filterName),  values = input.table[,ensembl.col],   mart = mart)       
                                                    
##head(results)



colnames(input.table)[ensembl.col] <- filterName
input.table <- merge(input.table, results, by=filterName, all.x=TRUE)

write.table(input.table, file=paste(file.name, "_symbol.csv", sep=""),  sep = "\t")



# 














