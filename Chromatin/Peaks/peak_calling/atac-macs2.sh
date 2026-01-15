#!/bin/sh
#$ -cwd
#$ -l h_vmem=4G
#$ -l h_rt=1:0:0
#$ -t 1-2

module load python/3.6.3
module load use.own
module load macs2

ip=$(sed -n "${SGE_TASK_ID}p" ip_list.txt)
exp=$(echo $ip | sed s/.bam//)

##PE data, narrow peaks
macs2 callpeak -t $ip -g mm -n $exp -q 0.05 -f BAMPE
