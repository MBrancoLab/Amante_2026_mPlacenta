#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G
#$ -t 1-39

module load python/3.6.3
module load use.own
module load localtools


te=$(sed -n "${SGE_TASK_ID}p" te_list.txt)

if echo $te | grep -q "L1"; then
	grep $te$'\t' ../mm_rmask_amante.bed | awk '$3-$2>5000' > $te.bed
fi

if echo $te | grep -q "IAP"; then
	grep $te$'\t' IAP_proviral_mm10.bed > $te.bed
fi

if echo $te | grep -q "RLTR4"; then
        grep $te$'\t' RLTR4_proviral_mm10.bed > $te.bed
fi


computeMatrix scale-regions -S H3K4me3*musculus*.bw IgG*musculus*.bw \
	-R $te.bed \
	-m 6000 \
	-a 2000 \
	-b 2000 \
	-bs 100 \
	-o $te.mat


plotProfile --matrixFile $te.mat \
        --outFileName $te.svg \
	--outFileNameData ${te}_trend.txt \
	--averageType "mean"
