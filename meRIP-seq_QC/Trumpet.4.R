library(Trumpet)
getwd()



## Collect the path of all the aligned MeRIP-seq data files in BAM format.
f1  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B3.IP.bam"
f2  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B5.IP.bam"
f3  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B6.IP.bam"
f4  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B7.IP.bam"
f5  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B8.IP.bam"
f6  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B9.IP.bam"
f7  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B12.IP.bam"
f8  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B13.IP.bam"
f9  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B16.IP.bam"
f10 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B20.IP.bam"
f11 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B26.IP.bam"
f12 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B28.IP.bam"
f13 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/IP/C_samples/5-finalBAM/3_STAR/B30.IP.bam"

if1  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B3.IN.bam"
if2  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B5.IN.bam"
if3  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B6.IN.bam"
if4  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B7.IN.bam"
if5  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B8.IN.bam"
if6  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B9.IN.bam"
if7  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B12.IN.bam"
if8  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B13.IN.bam"
if9  <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B16.IN.bam"
if10 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B20.IN.bam"
if11 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B26.IN.bam"
if12 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B28.IN.bam"
if13 <- "/media/yp/yongpeng16TB/PsychoENCODE/m6A-seq/Input/C_samples/5-finalBAM/3_STAR/B30.IN.bam"


ip_bam    <- c(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10)
input_bam <- c(if1,if2,if3,if4,if5,if6,if7,if8,if9,if10)

contrast_ip_bam    <- c(f11,f12,f13)
contrast_input_bam <- c(if11,if12,if13)

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



