#!/bin/sh
#$ -cwd
#$ -pe smp 4
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-5

bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
base=$(basename $bam | sed s/.bam//)

##merge gtf files
#~/stringtie-v2.2.3/stringtie --merge -o merged.gtf -G gencode.vM25.annotation.gtf Placenta_musculus-*.gtf 

##stringtie
mkdir ${base}_quant
~/stringtie-v2.2.3/stringtie --rf -e -B -G merged.gtf -o ${base}_quant/${base}.gtf -p 4 ${base}-sorted.bam
mv ${base}_quant/t_data.ctab ./${base}_quant.ctab
