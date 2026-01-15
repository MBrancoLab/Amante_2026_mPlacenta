#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-4


module load bedtools


file=$(sed -n "${SGE_TASK_ID}p" PL_bam_list.txt)
#file=$(sed -n "${SGE_TASK_ID}p" TS_bam_list.txt)
#file=$(sed -n "${SGE_TASK_ID}p" MP_bam_list.txt)
base_name=$(echo $file | sed s/.bam//)

##intersect
intersectBed -v -a ${base_name}.bam -b PL_blacklist.bed > ${base_name}.revised.bam
#intersectBed -v -a ${base_name}.bam -b TS_blacklist.bed > ${base_name}.revised.bam
#intersectBed -v -a ${base_name}.bam -b MP_blacklist.bed > ${base_name}.revised.bam
