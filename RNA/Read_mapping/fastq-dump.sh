#!/bin/sh
#$ -cwd
#$ -l h_vmem=2G
#$ -l h_rt=12:0:0
#$ -t 1-4

module load sratools/2.8.2-1

srr_url=$(sed -n "${SGE_TASK_ID}p" srr_list.txt)
srr=$(basename $srr_url)

##download
#wget $srr_url
fastq-dump --split-files --gzip $srr

##rename
mv ${srr}_1.fastq.gz E14_S1_L00${SGE_TASK_ID}_R1_001.fastq.gz
mv ${srr}_2.fastq.gz E14_S1_L00${SGE_TASK_ID}_R2_001.fastq.gz
mv ${srr}_3.fastq.gz E14_S1_L00${SGE_TASK_ID}_I1_001.fastq.gz

