#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_vmem=5G
#$ -j y
#$ -l h_rt=6:0:0
#$ -t 1-10

module load R
module load bedtools

peakFile=$(sed -n "${SGE_TASK_ID}p" peak_list.txt)

##for musculus
Rscript permutation_pipeline.R $peakFile mm_permutation_args.txt

##for pahari
#Rscript permutation_pipeline.R $peakFile mp_permutation_args.txt