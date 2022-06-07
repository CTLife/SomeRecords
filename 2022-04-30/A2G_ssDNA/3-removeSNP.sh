
mkdir  -p  3-removeSNP/GY6
bedtools intersect -v -a 2-sepAGCT/GY6.snp/Znormal_A.bed   -b mm39.vcf  >  3-removeSNP/GY6/A.bed
bedtools intersect -v -a 2-sepAGCT/GY6.snp/Znormal_G.bed   -b mm39.vcf  >  3-removeSNP/GY6/G.bed
bedtools intersect -v -a 2-sepAGCT/GY6.snp/Znormal_C.bed   -b mm39.vcf  >  3-removeSNP/GY6/C.bed
bedtools intersect -v -a 2-sepAGCT/GY6.snp/Znormal_T.bed   -b mm39.vcf  >  3-removeSNP/GY6/T.bed


mkdir  -p  3-removeSNP/MALBAC
bedtools intersect -v -a 2-sepAGCT/MALBAC.snp/Znormal_A.bed   -b mm39.vcf  >  3-removeSNP/MALBAC/A.bed
bedtools intersect -v -a 2-sepAGCT/MALBAC.snp/Znormal_G.bed   -b mm39.vcf  >  3-removeSNP/MALBAC/G.bed
bedtools intersect -v -a 2-sepAGCT/MALBAC.snp/Znormal_C.bed   -b mm39.vcf  >  3-removeSNP/MALBAC/C.bed
bedtools intersect -v -a 2-sepAGCT/MALBAC.snp/Znormal_T.bed   -b mm39.vcf  >  3-removeSNP/MALBAC/T.bed


mkdir  -p  3-removeSNP/MDA
bedtools intersect -v -a 2-sepAGCT/MDA.snp/Znormal_A.bed   -b mm39.vcf  >  3-removeSNP/MDA/A.bed
bedtools intersect -v -a 2-sepAGCT/MDA.snp/Znormal_G.bed   -b mm39.vcf  >  3-removeSNP/MDA/G.bed
bedtools intersect -v -a 2-sepAGCT/MDA.snp/Znormal_C.bed   -b mm39.vcf  >  3-removeSNP/MDA/C.bed
bedtools intersect -v -a 2-sepAGCT/MDA.snp/Znormal_T.bed   -b mm39.vcf  >  3-removeSNP/MDA/T.bed






