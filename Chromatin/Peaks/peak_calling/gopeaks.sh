#!/bin/sh
#$ -cwd
#$ -l h_vmem=64G
#$ -l h_rt=1:0:0
#$ -t 1-9

module load use.own
module load gopeaks
module load samtools

bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
exp=$(echo $bam | sed s/.bam//)


n=$(samtools view -c $bam)
cutoff=$(expr $n / 1000000)


##for narrow peaks
gopeaks -b $bam -o $exp -r $cutoff -t 200

##for broad peaks
#gopeaks -b $bam -o $exp -r $cutoff -t 2000 -l 200
