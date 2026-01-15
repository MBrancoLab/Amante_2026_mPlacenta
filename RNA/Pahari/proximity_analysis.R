library(tidyverse)
library(GenomicRanges)


##process read counts from Seqmonk

mm = read_tsv('../musculus_raw_counts.txt.gz') %>%
	select(Probe:`Probe Strand`,ID,Placenta_musculus_1:Placenta_musculus_3) %>%
	mutate(across(Placenta_musculus_1:Placenta_musculus_3, ~./sum(.)*1e6))

mp = read_tsv('../pahari_raw_counts.txt.gz') %>%
	select(Probe:`Probe Strand`,Placenta_pahari_2:Placenta_pahari_3) %>%
	mutate(across(Placenta_pahari_2:Placenta_pahari_3, ~./sum(.)*1e6))
colnames(mp)[1] = 'pahari_id'


##orthology info from biomart

ortho = read_tsv('mm.mp.orthologues.txt.gz')
ortho1 = filter(ortho, `Shrew mouse homology type`=='ortholog_one2one') %>%
	dplyr::select(`Gene stable ID`, `Shrew mouse gene stable ID`) %>%
	distinct()
colnames(ortho1) = c('ID', 'pahari_id')


##merge datasets
##remove sex chromosome data

mm.ortho = inner_join(mm, ortho1, by='ID') %>%
	relocate(pahari_id, .after=ID)

mm.mp = inner_join(mm.ortho, mp, by='pahari_id') %>%
	filter(Chromosome.x!='X', Chromosome.x!='Y', Chromosome.y!='X', Chromosome.y!='Y')


##get promoters

mm.rna = GRanges(seqnames=paste('chr',mm.mp$Chromosome.x,sep=''),
		ranges=IRanges(start=mm.mp$Start.x, end=mm.mp$End.x),
		strand=mm.mp$`Probe Strand.x`,
		id=mm.mp$ID,
		gene=mm.mp$Probe,
		musculus=rowMeans(log2(select(mm.mp, Placenta_musculus_1:Placenta_musculus_3)+0.05)),
		pahari=rowMeans(log2(select(mm.mp, Placenta_pahari_2:Placenta_pahari_3)+0.05)))

mp.rna = GRanges(seqnames=paste('chr',mm.mp$Chromosome.y,sep=''),
		ranges=IRanges(start=mm.mp$Start.y, end=mm.mp$End.y),
		strand=mm.mp$`Probe Strand.y`,
		id=mm.mp$pahari_id,
		gene=mm.mp$Probe,
		pahari=rowMeans(log2(select(mm.mp, Placenta_pahari_2:Placenta_pahari_3)+0.05)),
		musculus=rowMeans(log2(select(mm.mp, Placenta_musculus_1:Placenta_musculus_3)+0.05)))
				
mm.prom = promoters(mm.rna)
mp.prom = promoters(mp.rna)


##get nearest TE to TSS

nearestTE = function(tss.data, bed.file, breaks=c(0,1e4,5e4,1e5,1e9)) {
	bed = read.delim(bed.file, header=F)
	tes = GRanges(seqnames=bed$V1,
			ranges=IRanges(start=bed$V2, end=bed$V3),
			strand=bed$V6,
			id=paste(bed$V4,bed$V1,bed$V2,sep='_'))
	near = distanceToNearest(tss.data, tes, ignore.strand=TRUE)
	
	tss.data$nearest = tss.data$distance = rep(NA,length(tss.data))
	tss.data$nearest[near@from] = tes$id[near@to]
	tss.data$distance[near@from] = elementMetadata(near)$distance
	tss.data$dgroup = cut(tss.data$distance, breaks, include.lowest=T)
	
	return(tss.data)
}

iap = nearestTE(mm.prom, 'IAP_mm_K4me3.bed')
l1 = nearestTE(mm.prom, 'L1Md_mm_K4me3.bed')
epi = nearestTE(mm.prom, 'RMER10A_epiOrth_mm_K27ac.bed')
gen.mm = nearestTE(mm.prom, 'RMER10A_genOrth_mm_K27ac.bed')
gen.mp = nearestTE(mp.prom, 'RMER10A_genOrth_mp_K27ac.bed')


##filter

min.expr = -1

iap.df = filter(as.data.frame(iap), !is.na(distance), musculus>min.expr | pahari>min.expr)
l1.df = filter(as.data.frame(l1), !is.na(distance), musculus>min.expr | pahari>min.expr)
epi.df = filter(as.data.frame(epi), !is.na(distance), musculus>min.expr | pahari>min.expr)
gen.mm.df = filter(as.data.frame(gen.mm), !is.na(distance), musculus>min.expr | pahari>min.expr)
gen.mp.df = filter(as.data.frame(gen.mp), !is.na(distance), musculus>min.expr | pahari>min.expr)


##plot relative expression as a function of distance

plotFC = function(data) {
	p = ggplot(data, aes(x=subf, y=musculus-pahari, fill=dgroup)) +
		geom_boxplot(outlier.colour=NA) +
		theme_classic() +
		ylim(-3,3) +
		scale_fill_manual(values=colorRampPalette(c('lightblue','darkblue'))(4)) +
		geom_hline(yintercept=0, linetype='dashed') +
		xlab('') +
		ylab('Relative expression')
	
	print(p)
}	

quartz(w=4.5,h=3)
plotFC(rbind(add_column(iap.df, subf='IAP'), add_column(l1.df, subf='L1Md')))
plotFC(rbind(add_column(gen.mm.df, subf='mm'), add_column(gen.mp.df, subf='mp'), add_column(epi.df, subf='epi')))


##ANOVA

summary(aov(musculus-pahari ~ dgroup, data=iap.df))
summary(aov(musculus-pahari ~ dgroup, data=l1.df))
summary(aov(musculus-pahari ~ dgroup, data=epi.df))
summary(aov(musculus-pahari ~ dgroup, data=gen.mm.df))
summary(aov(musculus-pahari ~ dgroup, data=gen.mp.df))


##number of genes in each group

group_by(iap.df, dgroup) %>% summarise(n=length(start))
group_by(l1.df, dgroup) %>% summarise(n=length(start))
group_by(epi.df, dgroup) %>% summarise(n=length(start))
group_by(gen.mm.df, dgroup) %>% summarise(n=length(start))
group_by(gen.mp.df, dgroup) %>% summarise(n=length(start))


##write proximal genes to file

write_tsv(filter(iap.df, distance<1e4), 'IAP_proximal_genes.txt')
write_tsv(filter(gen.mm.df, distance<1e4), 'RMER10A_mm_proximal_genes.txt')
write_tsv(filter(gen.mp.df, distance<1e4), 'RMER10A_mp_proximal_genes.txt')
	

