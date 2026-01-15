#!/bin/sh
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=2:0:0
#$ -l h_vmem=8G
#$ -t 1-3

module load samtools
module load python
module load use.own
module load localtools


bam=$(sed -n "${SGE_TASK_ID}p" bam2bigwig_list.txt)
base=$(echo $bam | sed s/.bam//)


##sort and index bam
samtools sort -o $base-sorted.bam $bam
samtools index $base-sorted.bam


##make bigwig
bamCoverage --bam $base-sorted.bam -o $base.bw --binSize 200 --normalizeUsing CPM
