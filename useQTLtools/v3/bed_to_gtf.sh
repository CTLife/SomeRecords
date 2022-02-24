sort  -k1,1   -k2,2n   Peaks_27samples.bed     >  Peaks_27samples.sorted.bed 
bedToGenePred   Peaks_27samples.sorted.bed   Peaks_27samples.sorted.genePred
genePredToGtf   file   Peaks_27samples.sorted.genePred     Peaks_27samples.sorted.gtf

bgzip   Peaks_27samples.sorted.gtf
 



