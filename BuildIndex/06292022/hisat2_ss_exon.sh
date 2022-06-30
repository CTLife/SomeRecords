
prefix1=gencode.vM29.primary_assembly.annotation
hisat2_extract_splice_sites.py   $prefix1.gtf   >   $prefix1.ss
hisat2_extract_exons.py   $prefix1.gtf   >   $prefix1.exon

