#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=1:0:0
#$ -l h_vmem=2G

module load python/3.6.3
module load use.own
module load localtools



for te in RLTR13D5 RLTR13B3 RMER10A; do

computeMatrix reference-point -S H3K27ac_bTSCs_musculus.bw H3K27ac_GFPTSCs_musculus.bw H3K27ac_placenta_musculus.bw \
	-R $te.bed \
	--referencePoint "center" \
	-a 3000 \
	-b 3000 \
	-bs 50 \
	-o $te.mat

done

for te in RLTR13D5 RLTR13B3; do

plotHeatmap --matrixFile $te.mat \
	--outFileName $te.svg \
	--sortUsing max \
	--sortUsingSamples 1 \
	--heatmapHeight 12 \
	--heatmapWidth 3 \
	--whatToShow "heatmap and colorbar" \
	--colorMap "viridis"

done


te=RMER10A

plotHeatmap --matrixFile $te.mat \
        --outFileName $te.svg \
        --sortUsing max \
	--sortUsingSamples 3 \
        --heatmapHeight 12 \
        --heatmapWidth 3 \
        --whatToShow "heatmap and colorbar" \
	--colorMap "viridis" \
	--zMax 1.5
