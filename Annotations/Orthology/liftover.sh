loPath=~/Documents/liftOver
rmPath=~/Documents/Repeatmasker


te=L1MdA

grep $te ${rmPath}/mm_rmask_amante.bed > ${te}_Mm10.bed
grep $te ${rmPath}/mc_rmask_amante.bed > ${te}_GCF_900094665.1.bed
grep $te ${rmPath}/mp_rmask_amante.bed > ${te}_GCF_900095145.1.bed
grep $te ${rmPath}/rn_rmask_amante.bed > ${te}_Rn7.bed

echo "Mm10 annotated: $(wc -l ${te}_Mm10.bed)"
echo "Mm10 full-length: $(awk '$3-$2>5000' ${te}_Mm10.bed | wc -l)"

for sp in GCF_900094665.1 GCF_900095145.1 Rn7
do
	${loPath}/liftOver -minMatch=0.5 ${te}_Mm10.bed ${loPath}/mm10To${sp}.over.chain lifted.bed failed.bed
	intersectBed -wa -wb -a lifted.bed -b ${te}_${sp}.bed > intersect.bed

	#check that names match and takes only unique hits
	#we conservatively don't count regions that liftover but don't have a matching annotation
	awk '$4==$10' intersect.bed | cut -f 7-12 | uniq > orthologous.bed

	echo "${sp} annotated: $(wc -l ${te}_${sp}.bed)"
	echo "${sp} full-length: $(awk '$3-$2>5000' ${te}_${sp}.bed | wc -l)"
	echo "${sp} orthologous: $(wc -l orthologous.bed)"

	rm lifted.bed failed.bed intersect.bed orthologous.bed
done