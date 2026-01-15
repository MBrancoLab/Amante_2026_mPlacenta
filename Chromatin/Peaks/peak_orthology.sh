loPath=~/Documents/liftOver


for mark in H3K4me3 H3K4me1 H3K27ac
do
	#add ID column
	gzcat musculus_cuttag/${mark}_placenta.bed.gz| awk 'BEGIN{ FS = OFS = "\t" }{print $0, "mm_peak_"NR}' > mm.bed
	gzcat pahari_cuttag/${mark}_placenta_pahari.bed.gz| awk 'BEGIN{ FS = OFS = "\t" }{print $0, "mp_peak_"NR}' > mp.bed

	#musculus-specific regions
	${loPath}/liftOver -minMatch=0.5 mm.bed ${loPath}/mm10ToGCF_900095145.1.over.chain lifted.bed failed.bed
	grep -v \# failed.bed > ${mark}_nonOrthologous_mm.bed

	#epigenetically conserved peaks
	intersectBed -u -a lifted.bed -b mp.bed > ${mark}_epiOrthologous_mp.bed
	${loPath}/liftOver -minMatch=0.5 ${mark}_epiOrthologous_mp.bed ${loPath}/GCF_900095145.1ToMm10.over.chain ${mark}_epiOrthologous_mm.bed failed.bed

	#only genetically (but not epigenetically) conserved
	intersectBed -v -a lifted.bed -b mp.bed > ${mark}_genOrthologous_mp.bed
	${loPath}/liftOver -minMatch=0.5 ${mark}_genOrthologous_mp.bed ${loPath}/GCF_900095145.1ToMm10.over.chain ${mark}_genOrthologous_mm.bed failed.bed

	#pahari-specific regions
	${loPath}/liftOver -minMatch=0.5 mp.bed ${loPath}/GCF_900095145.1ToMm10.over.chain lifted.bed failed.bed
	grep -v \# failed.bed > ${mark}_nonOrthologous_mp.bed

	rm mm.bed mp.bed lifted.bed failed.bed
done