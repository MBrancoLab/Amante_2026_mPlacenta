#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=8G
#$ -t 1-5

module load use.own
module load hisat2
module load trimgalore
module load samtools

R1=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
R2=$(echo $R1 | sed s/_R1/_R2/)
base=$(echo $R1 | sed s/_R1.fq.gz//)

##trim
trim_galore --dont_gzip --paired $R1 $R2

##align
hisat2 --no-softclip --no-unal -x /data/Blizard-BrancoLab/Genomes/mm10/hisat/genome -1 ${base}_R1_val_1.fq -2 ${base}_R2_val_2.fq -S $base.sam
#hisat2 --no-softclip --no-unal -x /data/Blizard-BrancoLab/Genomes/pahari/hisat/pahari -1 ${base}_R1_val_1.fq -2 ${base}_R2_val_2.fq -S $base.sam

##filter based on MAPQ
samtools view -b -q 2 $base.sam > $base.bam
