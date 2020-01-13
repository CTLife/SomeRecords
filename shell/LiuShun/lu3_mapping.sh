#!/bin/bash

#SBATCH --job-name=m6Aseq
#SBATCH --output=lu3_mapping.log
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2000


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

#main setting
cutadapt_out_dir=/project2/mengjiechen/lulu/hek293/batch12172018
star_index_name_rep=rRNA
star_index_rep=/project2/chuanhe/shunliu/genomes/human/GRCh38/star_rRNA_index
star_index_name=GRCh38
star_index=/project2/chuanhe/shunliu/genomes/human/GRCh38/star_v27nERCC_index
star_out_dir=/project2/mengjiechen/shunliu/lulu/hek293_12172018/mapping
stringtie_out_dir=/project2/mengjiechen/shunliu/lulu/hek293_12172018/assembly_quant
gtf_file=/project2/chuanhe/shunliu/annotations/gencodev27/gencode.v27.annotation.gtf
bed12_file=/project2/chuanhe/shunliu/annotations/gencodev27/gencode.v27.annotation.bed12
chrom_size=/home/shunliu/software/bedtools-2.27.1/genomes/hg38.nopatch.chrom.sizes
tx_size=135000000
peak_out_dir=/project2/mengjiechen/shunliu/lulu/hek293_12172018/peak_calling
genome_fa=/project2/chuanhe/shunliu/genomes/human/GRCh38/chroms/genome.fa
ncpus=24


#control_input -- mapping, assembly and quantification
control_input=lu3_r1
control_input_read_length=50
control_input_star_strandness="F"
control_input_stringtie_strandness="--fr"

control_input_r1_fq_file=$cutadapt_out_dir/Lu-3.R1.clean.fa

echo -e "\n$control_input\n"

if [[ -d $star_out_dir/$control_input ]];then rm -rf $star_out_dir/$control_input;fi
mkdir -p $star_out_dir/$control_input

echo -e "\nstar align -- genome mapping\n"
time STAR --genomeDir $star_index --readFilesIn $control_input_r1_fq_file         \
    --outFileNamePrefix $star_out_dir/$control_input/ \
    --runThreadN $ncpus --genomeLoad NoSharedMemory     \
    --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1    \
    --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.06              \
    --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000         \
    --outSAMunmapped None --outFilterType BySJout --outSAMattributes NH HI AS NM MD \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 30000000000 \
    --outWigType bedGraph --outWigStrand Stranded --outSAMstrandField None --outWigReferencesPrefix chr

time samtools view -@ $ncpus -F 1548 -Shub $star_out_dir/$control_input/Aligned.sortedByCoord.out.bam | samtools sort -T $star_out_dir/$control_input -@ $ncpus -o $star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.bam -

echo -e "\nsamtools index\n"
time samtools index $star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.bam

echo -e "\nsamtools flagstat\n"
time samtools flagstat $star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.bam > $star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.flagstat.qc

echo -e "\nmark duplicates\n"
java -Xmx4G -jar $(which picard.jar) MarkDuplicates I=$star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.bam \
	O=$star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.dupmark.bam M=$star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.dupmark.dup.qc \
	VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

echo -e "\nPBC file output\n"
bamToBed -i $star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.dupmark.bam | awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$2,$3,$6}' | sort | uniq -c | \
	awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > \
	$star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.dupmark.pbc.qc
rm $star_out_dir/$control_input/$control_input.$star_index_name.align.sorted.dupmark.bam

echo -e "\nbedGraphToBigWig\n"
str[1]=+; str[2]=-;
for istr in 1 2
do
for imult in Unique UniqueMultiple
do
    #grep ^chr $star_out_dir/$control_input/Signal.$imult.str$istr.out.bg > $star_out_dir/$control_input/sig.tmp
    LC_COLLATE=C sort -k1,1 -k2,2n $star_out_dir/$control_input/Signal.$imult.str$istr.out.bg -o $star_out_dir/$control_input/Signal.$imult.str$istr.out.bg
    bedGraphToBigWig $star_out_dir/$control_input/Signal.$imult.str$istr.out.bg $chrom_size $star_out_dir/$control_input/Signal.$imult.strand${str[istr]}.bw
done
done

echo -e "\n$control_input mapping, assembly and quantification done!\n"

wait

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
