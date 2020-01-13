#!/bin/bash

#SBATCH --job-name=luseq
#SBATCH --output=Lu-1.clean.log
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=2000


######################
# Begin work section #
######################


echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)

raw_name=Lu-1_L3_A005
sample=Lu-1

time cutadapt -e 0.1 -n 1 -O 1 -q 10 -m 21 -a AGATCGGAAGAGCACACGTCT -A GATCGTCGGACTGTAGAACTC -o $sample.R1.adapter3trim.fastq.gz -p $sample.R2.adapter3trim.fastq.gz \
	$raw_name.R1.clean.fastq.gz $raw_name.R2.clean.fastq.gz > $sample.adapter3trim.log
time zcat $sample.R1.adapter3trim.fastq.gz | fastx_collapser -Q33 -i - -o $sample.R1.adapter3trim.collapse.fa
time zcat $sample.R2.adapter3trim.fastq.gz | fastx_collapser -Q33 -i - -o $sample.R2.adapter3trim.collapse.fa
time cutadapt -e 0.1 -n 1 -O 1 -m 16 -u 5 -u -11 -o $sample.R1.clean.fa $sample.R1.adapter3trim.collapse.fa > $sample.R1.barcoder5n3trim.log
time cutadapt -e 0.1 -n 1 -O 1 -m 16 -u 11 -u -5 -o $sample.R2.clean.fa $sample.R2.adapter3trim.collapse.fa > $sample.R2.barcoder5n3trim.log
rm $sample.R[12].adapter3trim.fastq.gz
rm $sample.R[12].adapter3trim.collapse.fa

echo -e "\nAll done\n"

echo Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
