library(tidyverse)
library(GenomicRanges)


##process read counts from Seqmonk

rpm = read_tsv('gene_raw_counts.txt') %>%
	mutate(across(KRAB_01:IAP_gRNA_03, ~./sum(.)*1e6))

ctrl = rowMeans(log2(rpm[,grep('KRAB',colnames(rpm))] + 0.05))
iap = rowMeans(log2(rpm[,grep('IAP_gRNA',colnames(rpm))] + 0.05))

rna = GRanges(seqnames=paste('chr',rpm$Chromosome,sep=''),
		ranges=IRanges(start=rpm$Start, end=rpm$End),
		strand=rpm$`Probe Strand`,
		id=rpm$ID,
		gene=rpm$Probe,
		ctrl, iap)
		
tss = promoters(rna)


##get nearest TE to TSS

nearestTE = function(tss.data, bed.file) {
	bed = read.delim(bed.file, header=F)
	tes = GRanges(seqnames=bed$V1,
			ranges=IRanges(start=bed$V2, end=bed$V3),
			strand=bed$V6,
			id=paste(bed$V4,bed$V1,bed$V2,sep='_'))
	near = distanceToNearest(tss.data, tes, ignore.strand=TRUE)
	
	tss.data$nearest = tss.data$distance = rep(NA,length(tss))
	tss.data$nearest[near@from] = tes$id[near@to]
	tss.data$distance[near@from] = elementMetadata(near)$distance
	
	return(tss.data)
}

ltr1 = nearestTE(tss, 'IAP_subsets/IAPLTR1_Mm.bed')
other = nearestTE(tss, 'IAP_subsets/notIAPLTR1_Mm.bed')

sgpos = nearestTE(tss, 'IAP_subsets/IAP_sgRNA_pos.bed')
sgneg = nearestTE(tss, 'IAP_subsets/IAP_sgRNA_neg.bed')


##filter

min.expr = -1
dmax = 1e4

ltr1.df = filter(as.data.frame(ltr1), distance<dmax, ctrl>min.expr | iap>min.expr)
other.df = filter(as.data.frame(other), distance<dmax, ctrl>min.expr | iap>min.expr)

sgpos.df = filter(as.data.frame(sgpos), distance<dmax, ctrl>min.expr | iap>min.expr)
sgneg.df = filter(as.data.frame(sgneg), distance<dmax, ctrl>min.expr | iap>min.expr)


##plot lfc for two groups

plotFC = function(df1, df2, name1, name2) {
	data = rbind(select(df1, gene, ctrl, iap) %>% add_column(subset=name1),
				select(df2, gene, ctrl, iap) %>% add_column(subset=name2))
	
	p = ggplot(data, aes(x=subset, y=iap-ctrl)) +
		geom_boxplot(outlier.colour=NA, fill='lightblue') +
		theme_classic() +
		ylim(-1.5,1.5) +
		geom_hline(yintercept=0, linetype='dashed') +
		xlab('') +
		ylab('Log2 fold change')
	
	print(p)
	return(data)
}

quartz(w=2, h=3)
ltr.data = plotFC(ltr1.df, other.df, 'IAPLTR1_Mm', 'Other IAPs')
sg.data = plotFC(sgpos.df, sgneg.df, 'sgRNA+', 'sgRNA-')


##stats

with(ltr.data, t.test(iap-ctrl ~ subset))
with(sg.data, t.test(iap-ctrl ~ subset))



