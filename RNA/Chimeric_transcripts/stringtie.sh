#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=4:0:0
#$ -l h_vmem=2G
#$ -t 1-5


bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
base=$(basename $bam | sed s/.bam//)

##sort
module load samtools
samtools sort -@ 4 -o ${base}-sorted.bam $bam

##stringtie
~/stringtie-v2.2.3/stringtie ${base}-sorted.bam --rf -G gencode.vM25.annotation.gtf -o $base.gtf -p 4

##get TSSs for transcripts that overlap gencode exons
module unload gcc
module load R
Rscript get_tss_known_exons.R $base.gtf ${base}-tss.bed

##intersect with Repeatmasker
module load bedtools
intersectBed -wa -wb -a ../mm_rmask_amante.bed -b ${base}-tss.bed > ${base}-chimTx.txt

