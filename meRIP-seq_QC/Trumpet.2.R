library(Trumpet)
getwd()



## Collect the path of all the aligned MeRIP-seq data files in BAM format.
f1 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/BA9_samples/5-finalBAM/3_STAR/B14.IP.bam"
f2 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/BA9_samples/5-finalBAM/3_STAR/B15.IP.bam"
f3 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/BA9_samples/5-finalBAM/3_STAR/B25.IP.bam"
f4 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/BA9_samples/5-finalBAM/3_STAR/B32.IP.bam"
f5 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/BA9_samples/5-finalBAM/3_STAR/B14.IN.bam"
f6 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/BA9_samples/5-finalBAM/3_STAR/B15.IN.bam"
f7 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/BA9_samples/5-finalBAM/3_STAR/B25.IN.bam"
f8 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/BA9_samples/5-finalBAM/3_STAR/B34.IP.bam"
f9 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/BA9_samples/5-finalBAM/3_STAR/B34.IN.bam"

ip_bam <- c(f1,f2,f3,f4)
input_bam <- c(f5,f6,f7)
contrast_ip_bam <- c(f8)
contrast_input_bam <- c(f9)

##We use GTF file in the following example.
gtf <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/hg38.ncbiRefSeq.gtf"


options(bitmapType="cairo")


##We can call the main function to generate the assessment report under current working directory.
trumpet_report <- Trumpet_report(IP_BAM = ip_bam, 
                                 Input_BAM = input_bam, 
                                 contrast_IP_BAM = contrast_ip_bam, 
                                 contrast_Input_BAM = contrast_input_bam, 
                                 condition1 = "untreated", 
                                 condition2 = "treat", 
                                 GENE_ANNO_GTF = gtf )



