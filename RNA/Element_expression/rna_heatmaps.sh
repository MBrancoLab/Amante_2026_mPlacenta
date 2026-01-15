#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools


for te in IAPLTR1_Mm L1MdTf_I RLTR4_Mm
do

if echo $te | grep -q "L1"; then
	grep $te$'\t' mm_rmask_amante.bed | awk '$3-$2>5000' > $te.bed
fi

if echo $te | grep -q "IAP"; then
	grep $te$'\t' IAP_proviral_mm10.bed > $te.bed
fi

if echo $te | grep -q "RLTR4"; then
        grep $te$'\t' RLTR4_proviral_mm10.bed > $te.bed
fi


computeMatrix scale-regions -S Placenta_musculus.bw GFPTSCs_musculus.bw bTSCs_musculus.bw \
	-R $te.bed \
	-m 6000 \
	-a 2000 \
	-b 2000 \
	-bs 100 \
	-o ${te}_rna.mat


plotHeatmap --matrixFile ${te}_rna.mat \
	--outFileName ${te}_rna.svg \
	--sortUsing sum \
	--sortUsingSamples 1 \
	--heatmapHeight 12 \
	--heatmapWidth 3 \
	--whatToShow "heatmap and colorbar" \
	--zMax 0.37 \
	--colorMap "viridis"

done
