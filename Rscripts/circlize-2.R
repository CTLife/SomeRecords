

library(circlize)

#################
DNAme.Hyper.NC_vs_IVFfresh = read.table(file="2_format/DMRs/1_NC_vs_IVFfresh/hyper.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hypo.NC_vs_IVFfresh  = read.table(file="2_format/DMRs/1_NC_vs_IVFfresh/hypo.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hyper.NC_vs_IVFfresh = cbind(DNAme.Hyper.NC_vs_IVFfresh[,1:3], "V4"=1) 
DNAme.Hypo.NC_vs_IVFfresh  = cbind(DNAme.Hypo.NC_vs_IVFfresh[,1:3], "V4"=-1) 
All.NC_vs_IVFfresh   = rbind(DNAme.Hyper.NC_vs_IVFfresh, DNAme.Hypo.NC_vs_IVFfresh)
dim(DNAme.Hyper.NC_vs_IVFfresh)
dim(DNAme.Hypo.NC_vs_IVFfresh)
dim(All.NC_vs_IVFfresh)
head(DNAme.Hyper.NC_vs_IVFfresh)
head(DNAme.Hypo.NC_vs_IVFfresh)

DNAme.Hyper.NC_vs_ICSIfresh = read.table(file="2_format/DMRs/2_NC_vs_ICSIfresh/hyper.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hypo.NC_vs_ICSIfresh  = read.table(file="2_format/DMRs/2_NC_vs_ICSIfresh/hypo.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hyper.NC_vs_ICSIfresh = cbind(DNAme.Hyper.NC_vs_ICSIfresh[,1:3], "V4"=1)
DNAme.Hypo.NC_vs_ICSIfresh  = cbind(DNAme.Hypo.NC_vs_ICSIfresh[,1:3], "V4"=-1)
All.NC_vs_ICSIfresh   = rbind(DNAme.Hyper.NC_vs_ICSIfresh, DNAme.Hypo.NC_vs_ICSIfresh)
dim(DNAme.Hyper.NC_vs_ICSIfresh)
dim(DNAme.Hypo.NC_vs_ICSIfresh)
dim(All.NC_vs_ICSIfresh)
head(DNAme.Hyper.NC_vs_ICSIfresh)
head(DNAme.Hypo.NC_vs_ICSIfresh)

DNAme.Hyper.NC_vs_IVFfrozen = read.table(file="2_format/DMRs/3_NC_vs_IVFfrozen/hyper.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hypo.NC_vs_IVFfrozen  = read.table(file="2_format/DMRs/3_NC_vs_IVFfrozen/hypo.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hyper.NC_vs_IVFfrozen = cbind(DNAme.Hyper.NC_vs_IVFfrozen[,1:3], "V4"=1)
DNAme.Hypo.NC_vs_IVFfrozen  = cbind(DNAme.Hypo.NC_vs_IVFfrozen[,1:3], "V4"=-1)
All.NC_vs_IVFfrozen   = rbind(DNAme.Hyper.NC_vs_IVFfrozen, DNAme.Hypo.NC_vs_IVFfrozen)
dim(DNAme.Hyper.NC_vs_IVFfrozen)
dim(DNAme.Hypo.NC_vs_IVFfrozen)
dim(All.NC_vs_IVFfrozen)
head(DNAme.Hyper.NC_vs_IVFfrozen)
head(DNAme.Hypo.NC_vs_IVFfrozen)

DNAme.Hyper.NC_vs_ICSIfrozen = read.table(file="2_format/DMRs/4_NC_vs_ICSIfrozen/hyper.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hypo.NC_vs_ICSIfrozen  = read.table(file="2_format/DMRs/4_NC_vs_ICSIfrozen/hypo.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
DNAme.Hyper.NC_vs_ICSIfrozen = cbind(DNAme.Hyper.NC_vs_ICSIfrozen[,1:3], "V4"=1)
DNAme.Hypo.NC_vs_ICSIfrozen  = cbind(DNAme.Hypo.NC_vs_ICSIfrozen[,1:3], "V4"=-1)
All.NC_vs_ICSIfrozen   = rbind(DNAme.Hyper.NC_vs_ICSIfrozen, DNAme.Hypo.NC_vs_ICSIfrozen)
dim(DNAme.Hyper.NC_vs_ICSIfrozen)
dim(DNAme.Hypo.NC_vs_ICSIfrozen)
dim(All.NC_vs_ICSIfrozen)
head(DNAme.Hyper.NC_vs_ICSIfrozen)
head(DNAme.Hypo.NC_vs_ICSIfrozen)




#############
H3K4me1.Hyper.NC_vs_IVFfresh = read.table(file="2_format/diffPeaks/H3K4me1/1_NC_vs_IVFfresh/IVFfresh-NC.Loss.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me1.Hypo.NC_vs_IVFfresh  = read.table(file="2_format/diffPeaks/H3K4me1/1_NC_vs_IVFfresh/IVFfresh-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me1.Hyper.NC_vs_IVFfresh = cbind(H3K4me1.Hyper.NC_vs_IVFfresh[,1:3], "V4"=1) 
H3K4me1.Hypo.NC_vs_IVFfresh  = cbind(H3K4me1.Hypo.NC_vs_IVFfresh[,1:3], "V4"=-1) 
All.NC_vs_IVFfresh   = rbind(H3K4me1.Hyper.NC_vs_IVFfresh, H3K4me1.Hypo.NC_vs_IVFfresh)
dim(H3K4me1.Hyper.NC_vs_IVFfresh)
dim(H3K4me1.Hypo.NC_vs_IVFfresh)
dim(All.NC_vs_IVFfresh)
head(H3K4me1.Hyper.NC_vs_IVFfresh)
head(H3K4me1.Hypo.NC_vs_IVFfresh)


H3K4me1.Hyper.NC_vs_IVFfrozen = read.table(file="2_format/diffPeaks/H3K4me1/3_NC_vs_IVFfrozen/IVFfrozen-NC.Loss.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me1.Hypo.NC_vs_IVFfrozen  = read.table(file="2_format/diffPeaks/H3K4me1/3_NC_vs_IVFfrozen/IVFfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me1.Hyper.NC_vs_IVFfrozen = cbind(H3K4me1.Hyper.NC_vs_IVFfrozen[,1:3], "V4"=1)
H3K4me1.Hypo.NC_vs_IVFfrozen  = cbind(H3K4me1.Hypo.NC_vs_IVFfrozen[,1:3], "V4"=-1)
All.NC_vs_IVFfrozen   = rbind(H3K4me1.Hyper.NC_vs_IVFfrozen, H3K4me1.Hypo.NC_vs_IVFfrozen)
dim(H3K4me1.Hyper.NC_vs_IVFfrozen)
dim(H3K4me1.Hypo.NC_vs_IVFfrozen)
dim(All.NC_vs_IVFfrozen)
head(H3K4me1.Hyper.NC_vs_IVFfrozen)
head(H3K4me1.Hypo.NC_vs_IVFfrozen)

H3K4me1.Hyper.NC_vs_ICSIfrozen = read.table(file="2_format/diffPeaks/H3K4me1/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Loss.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me1.Hypo.NC_vs_ICSIfrozen  = read.table(file="2_format/diffPeaks/H3K4me1/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me1.Hyper.NC_vs_ICSIfrozen = cbind(H3K4me1.Hyper.NC_vs_ICSIfrozen[,1:3], "V4"=1)
H3K4me1.Hypo.NC_vs_ICSIfrozen  = cbind(H3K4me1.Hypo.NC_vs_ICSIfrozen[,1:3], "V4"=-1)
All.NC_vs_ICSIfrozen   = rbind(H3K4me1.Hyper.NC_vs_ICSIfrozen, H3K4me1.Hypo.NC_vs_ICSIfrozen)
dim(H3K4me1.Hyper.NC_vs_ICSIfrozen)
dim(H3K4me1.Hypo.NC_vs_ICSIfrozen)
dim(All.NC_vs_ICSIfrozen)
head(H3K4me1.Hyper.NC_vs_ICSIfrozen)
head(H3K4me1.Hypo.NC_vs_ICSIfrozen)




#############
H3K4me3.Hyper.NC_vs_IVFfresh = read.table(file="2_format/diffPeaks/H3K4me3/1_NC_vs_IVFfresh/IVFfresh-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hypo.NC_vs_IVFfresh  = read.table(file="2_format/diffPeaks/H3K4me3/1_NC_vs_IVFfresh/IVFfresh-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hyper.NC_vs_IVFfresh = cbind(H3K4me3.Hyper.NC_vs_IVFfresh[,1:3], "V4"=1) 
H3K4me3.Hypo.NC_vs_IVFfresh  = cbind(H3K4me3.Hypo.NC_vs_IVFfresh[,1:3], "V4"=-1) 
All.NC_vs_IVFfresh   = rbind(H3K4me3.Hyper.NC_vs_IVFfresh, H3K4me3.Hypo.NC_vs_IVFfresh)
dim(H3K4me3.Hyper.NC_vs_IVFfresh)
dim(H3K4me3.Hypo.NC_vs_IVFfresh)
dim(All.NC_vs_IVFfresh)
head(H3K4me3.Hyper.NC_vs_IVFfresh)
head(H3K4me3.Hypo.NC_vs_IVFfresh)

H3K4me3.Hyper.NC_vs_ICSIfresh = read.table(file="2_format/diffPeaks/H3K4me3/2_NC_vs_ICSIfresh/ICSIfresh-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hypo.NC_vs_ICSIfresh  = read.table(file="2_format/diffPeaks/H3K4me3/2_NC_vs_ICSIfresh/ICSIfresh-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hyper.NC_vs_ICSIfresh = cbind(H3K4me3.Hyper.NC_vs_ICSIfresh[,1:3], "V4"=1)
H3K4me3.Hypo.NC_vs_ICSIfresh  = cbind(H3K4me3.Hypo.NC_vs_ICSIfresh[,1:3], "V4"=-1)
All.NC_vs_ICSIfresh   = rbind(H3K4me3.Hyper.NC_vs_ICSIfresh, H3K4me3.Hypo.NC_vs_ICSIfresh)
dim(H3K4me3.Hyper.NC_vs_ICSIfresh)
dim(H3K4me3.Hypo.NC_vs_ICSIfresh)
dim(All.NC_vs_ICSIfresh)
head(H3K4me3.Hyper.NC_vs_ICSIfresh)
head(H3K4me3.Hypo.NC_vs_ICSIfresh)

H3K4me3.Hyper.NC_vs_IVFfrozen = read.table(file="2_format/diffPeaks/H3K4me3/3_NC_vs_IVFfrozen/IVFfrozen-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hypo.NC_vs_IVFfrozen  = read.table(file="2_format/diffPeaks/H3K4me3/3_NC_vs_IVFfrozen/IVFfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hyper.NC_vs_IVFfrozen = cbind(H3K4me3.Hyper.NC_vs_IVFfrozen[,1:3], "V4"=1)
H3K4me3.Hypo.NC_vs_IVFfrozen  = cbind(H3K4me3.Hypo.NC_vs_IVFfrozen[,1:3], "V4"=-1)
All.NC_vs_IVFfrozen   = rbind(H3K4me3.Hyper.NC_vs_IVFfrozen, H3K4me3.Hypo.NC_vs_IVFfrozen)
dim(H3K4me3.Hyper.NC_vs_IVFfrozen)
dim(H3K4me3.Hypo.NC_vs_IVFfrozen)
dim(All.NC_vs_IVFfrozen)
head(H3K4me3.Hyper.NC_vs_IVFfrozen)
head(H3K4me3.Hypo.NC_vs_IVFfrozen)

H3K4me3.Hyper.NC_vs_ICSIfrozen = read.table(file="2_format/diffPeaks/H3K4me3/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hypo.NC_vs_ICSIfrozen  = read.table(file="2_format/diffPeaks/H3K4me3/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K4me3.Hyper.NC_vs_ICSIfrozen = cbind(H3K4me3.Hyper.NC_vs_ICSIfrozen[,1:3], "V4"=1)
H3K4me3.Hypo.NC_vs_ICSIfrozen  = cbind(H3K4me3.Hypo.NC_vs_ICSIfrozen[,1:3], "V4"=-1)
All.NC_vs_ICSIfrozen   = rbind(H3K4me3.Hyper.NC_vs_ICSIfrozen, H3K4me3.Hypo.NC_vs_ICSIfrozen)
dim(H3K4me3.Hyper.NC_vs_ICSIfrozen)
dim(H3K4me3.Hypo.NC_vs_ICSIfrozen)
dim(All.NC_vs_ICSIfrozen)
head(H3K4me3.Hyper.NC_vs_ICSIfrozen)
head(H3K4me3.Hypo.NC_vs_ICSIfrozen)



#############
H3K27me3.Hyper.NC_vs_IVFfresh = read.table(file="2_format/diffPeaks/H3K27me3/1_NC_vs_IVFfresh/IVFfresh-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hypo.NC_vs_IVFfresh  = read.table(file="2_format/diffPeaks/H3K27me3/1_NC_vs_IVFfresh/IVFfresh-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hyper.NC_vs_IVFfresh = cbind(H3K27me3.Hyper.NC_vs_IVFfresh[,1:3], "V4"=1) 
H3K27me3.Hypo.NC_vs_IVFfresh  = cbind(H3K27me3.Hypo.NC_vs_IVFfresh[,1:3], "V4"=-1) 
All.NC_vs_IVFfresh   = rbind(H3K27me3.Hyper.NC_vs_IVFfresh, H3K27me3.Hypo.NC_vs_IVFfresh)
dim(H3K27me3.Hyper.NC_vs_IVFfresh)
dim(H3K27me3.Hypo.NC_vs_IVFfresh)
dim(All.NC_vs_IVFfresh)
head(H3K27me3.Hyper.NC_vs_IVFfresh)
head(H3K27me3.Hypo.NC_vs_IVFfresh)

H3K27me3.Hyper.NC_vs_ICSIfresh = read.table(file="2_format/diffPeaks/H3K27me3/2_NC_vs_ICSIfresh/ICSIfresh-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hypo.NC_vs_ICSIfresh  = read.table(file="2_format/diffPeaks/H3K27me3/2_NC_vs_ICSIfresh/ICSIfresh-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hyper.NC_vs_ICSIfresh = cbind(H3K27me3.Hyper.NC_vs_ICSIfresh[,1:3], "V4"=1)
H3K27me3.Hypo.NC_vs_ICSIfresh  = cbind(H3K27me3.Hypo.NC_vs_ICSIfresh[,1:3], "V4"=-1)
All.NC_vs_ICSIfresh   = rbind(H3K27me3.Hyper.NC_vs_ICSIfresh, H3K27me3.Hypo.NC_vs_ICSIfresh)
dim(H3K27me3.Hyper.NC_vs_ICSIfresh)
dim(H3K27me3.Hypo.NC_vs_ICSIfresh)
dim(All.NC_vs_ICSIfresh)
head(H3K27me3.Hyper.NC_vs_ICSIfresh)
head(H3K27me3.Hypo.NC_vs_ICSIfresh)

H3K27me3.Hyper.NC_vs_IVFfrozen = read.table(file="2_format/diffPeaks/H3K27me3/3_NC_vs_IVFfrozen/IVFfrozen-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hypo.NC_vs_IVFfrozen  = read.table(file="2_format/diffPeaks/H3K27me3/3_NC_vs_IVFfrozen/IVFfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hyper.NC_vs_IVFfrozen = cbind(H3K27me3.Hyper.NC_vs_IVFfrozen[,1:3], "V4"=1)
H3K27me3.Hypo.NC_vs_IVFfrozen  = cbind(H3K27me3.Hypo.NC_vs_IVFfrozen[,1:3], "V4"=-1)
All.NC_vs_IVFfrozen   = rbind(H3K27me3.Hyper.NC_vs_IVFfrozen, H3K27me3.Hypo.NC_vs_IVFfrozen)
dim(H3K27me3.Hyper.NC_vs_IVFfrozen)
dim(H3K27me3.Hypo.NC_vs_IVFfrozen)
dim(All.NC_vs_IVFfrozen)
head(H3K27me3.Hyper.NC_vs_IVFfrozen)
head(H3K27me3.Hypo.NC_vs_IVFfrozen)

H3K27me3.Hyper.NC_vs_ICSIfrozen = read.table(file="2_format/diffPeaks/H3K27me3/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hypo.NC_vs_ICSIfrozen  = read.table(file="2_format/diffPeaks/H3K27me3/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27me3.Hyper.NC_vs_ICSIfrozen = cbind(H3K27me3.Hyper.NC_vs_ICSIfrozen[,1:3], "V4"=1)
H3K27me3.Hypo.NC_vs_ICSIfrozen  = cbind(H3K27me3.Hypo.NC_vs_ICSIfrozen[,1:3], "V4"=-1)
All.NC_vs_ICSIfrozen   = rbind(H3K27me3.Hyper.NC_vs_ICSIfrozen, H3K27me3.Hypo.NC_vs_ICSIfrozen)
dim(H3K27me3.Hyper.NC_vs_ICSIfrozen)
dim(H3K27me3.Hypo.NC_vs_ICSIfrozen)
dim(All.NC_vs_ICSIfrozen)
head(H3K27me3.Hyper.NC_vs_ICSIfrozen)
head(H3K27me3.Hypo.NC_vs_ICSIfrozen)




#############
H3K27ac.Hyper.NC_vs_IVFfresh = read.table(file="2_format/diffPeaks/H3K27ac/1_NC_vs_IVFfresh/IVFfresh-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27ac.Hypo.NC_vs_IVFfresh  = read.table(file="2_format/diffPeaks/H3K27ac/1_NC_vs_IVFfresh/IVFfresh-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27ac.Hyper.NC_vs_IVFfresh = cbind(H3K27ac.Hyper.NC_vs_IVFfresh[,1:3], "V4"=1) 
H3K27ac.Hypo.NC_vs_IVFfresh  = cbind(H3K27ac.Hypo.NC_vs_IVFfresh[,1:3], "V4"=-1) 
All.NC_vs_IVFfresh   = rbind(H3K27ac.Hyper.NC_vs_IVFfresh, H3K27ac.Hypo.NC_vs_IVFfresh)
dim(H3K27ac.Hyper.NC_vs_IVFfresh)
dim(H3K27ac.Hypo.NC_vs_IVFfresh)
dim(All.NC_vs_IVFfresh)
head(H3K27ac.Hyper.NC_vs_IVFfresh)
head(H3K27ac.Hypo.NC_vs_IVFfresh)

H3K27ac.Hyper.NC_vs_IVFfrozen = read.table(file="2_format/diffPeaks/H3K27ac/3_NC_vs_IVFfrozen/IVFfrozen-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27ac.Hypo.NC_vs_IVFfrozen  = read.table(file="2_format/diffPeaks/H3K27ac/3_NC_vs_IVFfrozen/IVFfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27ac.Hyper.NC_vs_IVFfrozen = cbind(H3K27ac.Hyper.NC_vs_IVFfrozen[,1:3], "V4"=1)
H3K27ac.Hypo.NC_vs_IVFfrozen  = cbind(H3K27ac.Hypo.NC_vs_IVFfrozen[,1:3], "V4"=-1)
All.NC_vs_IVFfrozen   = rbind(H3K27ac.Hyper.NC_vs_IVFfrozen, H3K27ac.Hypo.NC_vs_IVFfrozen)
dim(H3K27ac.Hyper.NC_vs_IVFfrozen)
dim(H3K27ac.Hypo.NC_vs_IVFfrozen)
dim(All.NC_vs_IVFfrozen)
head(H3K27ac.Hyper.NC_vs_IVFfrozen)
head(H3K27ac.Hypo.NC_vs_IVFfrozen)

H3K27ac.Hyper.NC_vs_ICSIfrozen = read.table(file="2_format/diffPeaks/H3K27ac/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Loss.bed", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27ac.Hypo.NC_vs_ICSIfrozen  = read.table(file="2_format/diffPeaks/H3K27ac/4_NC_vs_ICSIfrozen/ICSIfrozen-NC.Gain.bed",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
H3K27ac.Hyper.NC_vs_ICSIfrozen = cbind(H3K27ac.Hyper.NC_vs_ICSIfrozen[,1:3], "V4"=1)
H3K27ac.Hypo.NC_vs_ICSIfrozen  = cbind(H3K27ac.Hypo.NC_vs_ICSIfrozen[,1:3], "V4"=-1)
All.NC_vs_ICSIfrozen   = rbind(H3K27ac.Hyper.NC_vs_ICSIfrozen, H3K27ac.Hypo.NC_vs_ICSIfrozen)
dim(H3K27ac.Hyper.NC_vs_ICSIfrozen)
dim(H3K27ac.Hypo.NC_vs_ICSIfrozen)
dim(All.NC_vs_ICSIfrozen)
head(H3K27ac.Hyper.NC_vs_ICSIfrozen)
head(H3K27ac.Hypo.NC_vs_ICSIfrozen)





#################
RNA.Hyper.NC_vs_IVFfresh = read.table(file="2_format/DEGs2/1_NC_vs_IVFfresh/up.txt", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hypo.NC_vs_IVFfresh  = read.table(file="2_format/DEGs2/1_NC_vs_IVFfresh/down.txt",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hyper.NC_vs_IVFfresh = cbind(RNA.Hyper.NC_vs_IVFfresh[,1:3], "V4"=1) 
RNA.Hypo.NC_vs_IVFfresh  = cbind(RNA.Hypo.NC_vs_IVFfresh[,1:3], "V4"=-1) 
All.NC_vs_IVFfresh   = rbind(RNA.Hyper.NC_vs_IVFfresh, RNA.Hypo.NC_vs_IVFfresh)
dim(RNA.Hyper.NC_vs_IVFfresh)
dim(RNA.Hypo.NC_vs_IVFfresh)
dim(All.NC_vs_IVFfresh)
head(RNA.Hyper.NC_vs_IVFfresh)
head(RNA.Hypo.NC_vs_IVFfresh)

RNA.Hyper.NC_vs_ICSIfresh = read.table(file="2_format/DEGs2/2_NC_vs_ICSIfresh/up.txt", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hypo.NC_vs_ICSIfresh  = read.table(file="2_format/DEGs2/2_NC_vs_ICSIfresh/down.txt",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hyper.NC_vs_ICSIfresh = cbind(RNA.Hyper.NC_vs_ICSIfresh[,1:3], "V4"=1)
RNA.Hypo.NC_vs_ICSIfresh  = cbind(RNA.Hypo.NC_vs_ICSIfresh[,1:3], "V4"=-1)
All.NC_vs_ICSIfresh   = rbind(RNA.Hyper.NC_vs_ICSIfresh, RNA.Hypo.NC_vs_ICSIfresh)
dim(RNA.Hyper.NC_vs_ICSIfresh)
dim(RNA.Hypo.NC_vs_ICSIfresh)
dim(All.NC_vs_ICSIfresh)
head(RNA.Hyper.NC_vs_ICSIfresh)
head(RNA.Hypo.NC_vs_ICSIfresh)

RNA.Hyper.NC_vs_IVFfrozen = read.table(file="2_format/DEGs2/3_NC_vs_IVFfrozen/up.txt", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hypo.NC_vs_IVFfrozen  = read.table(file="2_format/DEGs2/3_NC_vs_IVFfrozen/down.txt",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hyper.NC_vs_IVFfrozen = cbind(RNA.Hyper.NC_vs_IVFfrozen[,1:3], "V4"=1)
RNA.Hypo.NC_vs_IVFfrozen  = cbind(RNA.Hypo.NC_vs_IVFfrozen[,1:3], "V4"=-1)
All.NC_vs_IVFfrozen   = rbind(RNA.Hyper.NC_vs_IVFfrozen, RNA.Hypo.NC_vs_IVFfrozen)
dim(RNA.Hyper.NC_vs_IVFfrozen)
dim(RNA.Hypo.NC_vs_IVFfrozen)
dim(All.NC_vs_IVFfrozen)
head(RNA.Hyper.NC_vs_IVFfrozen)
head(RNA.Hypo.NC_vs_IVFfrozen)

RNA.Hyper.NC_vs_ICSIfrozen = read.table(file="2_format/DEGs2/4_NC_vs_ICSIfrozen/up.txt", header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hypo.NC_vs_ICSIfrozen  = read.table(file="2_format/DEGs2/4_NC_vs_ICSIfrozen/down.txt",  header = FALSE, sep = "\t", quote = "",  comment.char = "#")
RNA.Hyper.NC_vs_ICSIfrozen = cbind(RNA.Hyper.NC_vs_ICSIfrozen[,1:3], "V4"=1)
RNA.Hypo.NC_vs_ICSIfrozen  = cbind(RNA.Hypo.NC_vs_ICSIfrozen[,1:3], "V4"=-1)
All.NC_vs_ICSIfrozen   = rbind(RNA.Hyper.NC_vs_ICSIfrozen, RNA.Hypo.NC_vs_ICSIfrozen)
dim(RNA.Hyper.NC_vs_ICSIfrozen)
dim(RNA.Hypo.NC_vs_ICSIfrozen)
dim(All.NC_vs_ICSIfrozen)
head(RNA.Hyper.NC_vs_ICSIfrozen)
head(RNA.Hypo.NC_vs_ICSIfrozen)






myChrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
           "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
           "chr17", "chr18", "chr19", "chr20", "chr21", "chr22" )


Part1_g    = "Results"
if( ! file.exists(Part1_g)    ) { dir.create(Part1_g,    recursive = TRUE) }


pdf(file = paste("Results", "1.pdf", sep="/"), width=20, height=20 )
circos.par(start.degree=90, gap.degree=c(rep(0.2, 21), 10),  track.margin=c(0, 0.01),  track.height= 0.02, cell.padding=c(0, 0, 0, 0) )
circos.initializeWithIdeogram(species = "hg38", chromosome.index=myChrs  )
text(0, 0, NULL, cex = 5)

bed_list1 = list(RNA.Hyper.NC_vs_IVFfresh, RNA.Hypo.NC_vs_IVFfresh)
circos.genomicRainfall(bed_list1, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0) )
bed_list2 = list(RNA.Hyper.NC_vs_IVFfrozen, RNA.Hypo.NC_vs_IVFfrozen)
circos.genomicRainfall(bed_list2, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list3 = list(RNA.Hyper.NC_vs_ICSIfresh, RNA.Hypo.NC_vs_ICSIfresh)
circos.genomicRainfall(bed_list3, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list4 = list(RNA.Hyper.NC_vs_ICSIfrozen, RNA.Hypo.NC_vs_ICSIfrozen)
circos.genomicRainfall(bed_list4, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0.01) )
circos.par(track.height= 0.05 )

bed_list1 = list(DNAme.Hyper.NC_vs_IVFfresh, DNAme.Hypo.NC_vs_IVFfresh)
circos.genomicRainfall(bed_list1, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 8) )
circos.par(track.margin=c(0, 0) )
bed_list2 = list(DNAme.Hyper.NC_vs_IVFfrozen, DNAme.Hypo.NC_vs_IVFfrozen)
circos.genomicRainfall(bed_list2, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 8) )
bed_list3 = list(DNAme.Hyper.NC_vs_ICSIfresh, DNAme.Hypo.NC_vs_ICSIfresh)
circos.genomicRainfall(bed_list3, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 8) )
bed_list4 = list(DNAme.Hyper.NC_vs_ICSIfrozen, DNAme.Hypo.NC_vs_ICSIfrozen)
circos.genomicRainfall(bed_list4, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 8) )
circos.par(track.margin=c(0, 0.01) )

bed_list1 = list(H3K4me3.Hyper.NC_vs_IVFfresh, H3K4me3.Hypo.NC_vs_IVFfresh)
circos.genomicRainfall(bed_list1, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0) )
bed_list2 = list(H3K4me3.Hyper.NC_vs_IVFfrozen, H3K4me3.Hypo.NC_vs_IVFfrozen)
circos.genomicRainfall(bed_list2, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list3 = list(H3K4me3.Hyper.NC_vs_ICSIfresh, H3K4me3.Hypo.NC_vs_ICSIfresh)
circos.genomicRainfall(bed_list3, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list4 = list(H3K4me3.Hyper.NC_vs_ICSIfrozen, H3K4me3.Hypo.NC_vs_ICSIfrozen)
circos.genomicRainfall(bed_list4, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0.01) )
circos.par(track.height= 0.02 )

bed_list1 = list(H3K27ac.Hyper.NC_vs_IVFfresh, H3K27ac.Hypo.NC_vs_IVFfresh)
circos.genomicRainfall(bed_list1, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0) )
bed_list2 = list(H3K27ac.Hyper.NC_vs_IVFfrozen, H3K27ac.Hypo.NC_vs_IVFfrozen)
circos.genomicRainfall(bed_list2, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list3 = list()
circos.genomicRainfall(bed_list3, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list4 = list(H3K27ac.Hyper.NC_vs_ICSIfrozen, H3K27ac.Hypo.NC_vs_ICSIfrozen)
circos.genomicRainfall(bed_list4, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0.01) )

bed_list1 = list(H3K4me1.Hyper.NC_vs_IVFfresh, H3K4me1.Hypo.NC_vs_IVFfresh)
circos.genomicRainfall(bed_list1, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0) )
bed_list2 = list(H3K4me1.Hyper.NC_vs_IVFfrozen, H3K4me1.Hypo.NC_vs_IVFfrozen)
circos.genomicRainfall(bed_list2, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list3 = list()
circos.genomicRainfall(bed_list3, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list4 = list(H3K4me1.Hyper.NC_vs_ICSIfrozen, H3K4me1.Hypo.NC_vs_ICSIfrozen)
circos.genomicRainfall(bed_list4, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0.01) )

bed_list1 = list(H3K27me3.Hyper.NC_vs_IVFfresh, H3K27me3.Hypo.NC_vs_IVFfresh)
circos.genomicRainfall(bed_list1, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0) )
bed_list2 = list(H3K27me3.Hyper.NC_vs_IVFfrozen, H3K27me3.Hypo.NC_vs_IVFfrozen)
circos.genomicRainfall(bed_list2, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list3 = list(H3K27me3.Hyper.NC_vs_ICSIfresh, H3K27me3.Hypo.NC_vs_ICSIfresh)
circos.genomicRainfall(bed_list3, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
bed_list4 = list(H3K27me3.Hyper.NC_vs_ICSIfrozen, H3K27me3.Hypo.NC_vs_ICSIfrozen)
circos.genomicRainfall(bed_list4, pch = 16, cex = 0.4, col = c("red4", "cyan4"),  ylim = c(0, 9) )
circos.par(track.margin=c(0, 0.01) )

circos.clear()
dev.off()





