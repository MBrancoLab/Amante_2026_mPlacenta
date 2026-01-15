#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-8


module load bowtie2
module load trimgalore
module load samtools


file=$(sed -n "${SGE_TASK_ID}p" fastq_list.txt)
base_name=$(echo $file | sed s/_1.fq.gz//)

##trim
trim_galore --dont_gzip --paired ${base_name}_1.fq.gz ${base_name}_2.fq.gz

##map
bowtie2 -p 4 --no-unal -x /data/Blizard-BrancoLab/Genomes/mm10/mm10 -1 ${base_name}_1_val_1.fq -2 ${base_name}_2_val_2.fq -S ${base_name}.sam
#bowtie2 -p 4 --no-unal -x /data/Blizard-BrancoLab/Genomes/pahari/pahari -1 ${base_name}_1_val_1.fq -2 ${base_name}_2_val_2.fq -S ${base_name}.sam

##filter based on MAPQ
samtools view -b -q 2 $base_name.sam > $base_name.bam