library(Trumpet)
getwd()


## Collect the path of all the aligned MeRIP-seq data files in BAM format.
f1 <- system.file("extdata", "IP1.bam", package="Trumpet")
f2 <- system.file("extdata", "IP2.bam", package="Trumpet")
f3 <- system.file("extdata", "IP3.bam", package="Trumpet")
f4 <- system.file("extdata", "IP4.bam", package="Trumpet")
f5 <- system.file("extdata", "Input1.bam", package="Trumpet")
f6 <- system.file("extdata", "Input2.bam", package="Trumpet")
f7 <- system.file("extdata", "Input3.bam", package="Trumpet")
f8 <- system.file("extdata", "treated_IP1.bam", package="Trumpet")
f9 <- system.file("extdata", "treated_Input1.bam", package="Trumpet")

ip_bam <- c(f1,f2,f3,f4)
input_bam <- c(f5,f6,f7)
contrast_ip_bam <- c(f8)
contrast_input_bam <- c(f9)

##We use GTF file in the following example.
gtf <- system.file("extdata", "hg19toy.gtf", package="Trumpet")

options(bitmapType="cairo")



##We can call the main function to generate the assessment report under current working directory.
trumpet_report <- Trumpet_report(IP_BAM = ip_bam, 
                                 Input_BAM = input_bam, 
                                 contrast_IP_BAM = contrast_ip_bam, 
                                 contrast_Input_BAM = contrast_input_bam, 
                                 condition1 = "untreated", 
                                 condition2 = "treat", 
                                 GENE_ANNO_GTF = gtf 
                                 )


                                 
                                 
                                 
                                 
                                 
