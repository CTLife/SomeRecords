out="2-fasta"
mkdir -p $out

bedtools getfasta  -fo $out/TP.3UTR.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/TP.3UTR.bed                


bedtools getfasta  -fo $out/TP.5UTR.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/TP.5UTR.bed                



bedtools getfasta  -fo $out/TP.Exon.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/TP.Exon.bed                


bedtools getfasta  -fo $out/TP.all.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/TP.all.bed                




bedtools getfasta  -fo $out/FP.3UTR.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/FP.3UTR.bed                


bedtools getfasta  -fo $out/FP.5UTR.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/FP.5UTR.bed                



bedtools getfasta  -fo $out/FP.Exon.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/FP.Exon.bed                


bedtools getfasta  -fo $out/FP.all.fasta    -name  -s    -fi hg38.fa   \
                   -fullHeader   -bed 1-bed/FP.all.bed                





