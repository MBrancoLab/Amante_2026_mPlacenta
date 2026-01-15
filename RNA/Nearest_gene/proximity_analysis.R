library(tidyverse)
library(GenomicRanges)


##process read counts from Seqmonk

rpm = read_tsv('../musculus_raw_counts.txt.gz') %>%
	mutate(across(bTSCs_musculus:Placenta_musculus_3, ~./sum(.)*1e6))

tsc = log2(rpm$GFPTSCs_musculus + 0.05)
plac = rowMeans(log2(rpm[,grep('Placenta',colnames(rpm))] + 0.05))
plac = plac + quantile(tsc,0.75) - quantile(plac,0.75)  ##adjust skewed distributions

rna = GRanges(seqnames=paste('chr',rpm$Chromosome,sep=''),
		ranges=IRanges(start=rpm$Start, end=rpm$End),
		strand=rpm$`Probe Strand`,
		id=rpm$ID,
		gene=rpm$Probe,
		tsc, plac)
		
tss = promoters(rna)


##get nearest TE to TSS

nearestTE = function(tss.data, bed.file, breaks=c(0,1e4,5e4,1e5,1e9)) {
	bed = read.delim(bed.file, header=F)
	tes = GRanges(seqnames=bed$V1,
			ranges=IRanges(start=bed$V2, end=bed$V3),
			strand=bed$V6,
			id=paste(bed$V4,bed$V1,bed$V2,sep='_'))
	near = distanceToNearest(tss.data, tes, ignore.strand=TRUE)
	
	tss.data$nearest = tss.data$distance = rep(NA,length(tss))
	tss.data$nearest[near@from] = tes$id[near@to]
	tss.data$distance[near@from] = elementMetadata(near)$distance
	tss.data$dgroup = cut(tss.data$distance, breaks, include.lowest=T)
	
	return(tss.data)
}

rltr = nearestTE(tss, 'TSC_RLTR13D5.bed')
rmer = nearestTE(tss, 'Placenta_RMER10A.bed')
iap.c = nearestTE(tss, 'Common_IAP.bed')
iap.p = nearestTE(tss, 'Placenta_IAP.bed')
l1.c = nearestTE(tss, 'Common_L1Md.bed')
l1.p = nearestTE(tss, 'Placenta_L1Md.bed')


##filter

min.expr = -1

rltr.df = filter(as.data.frame(rltr), !is.na(distance), tsc>min.expr | plac>min.expr)
rmer.df = filter(as.data.frame(rmer), !is.na(distance), tsc>min.expr | plac>min.expr)
iapc.df = filter(as.data.frame(iap.c), !is.na(distance), tsc>min.expr | plac>min.expr)
iapp.df = filter(as.data.frame(iap.p), !is.na(distance), tsc>min.expr | plac>min.expr)
l1c.df = filter(as.data.frame(l1.c), !is.na(distance), tsc>min.expr | plac>min.expr)
l1p.df = filter(as.data.frame(l1.p), !is.na(distance), tsc>min.expr | plac>min.expr)


##plot relative expression as a function of distance

plotFC = function(data) {
	p = ggplot(data, aes(x=subf, y=plac-tsc, fill=dgroup)) +
		geom_boxplot(outlier.colour=NA) +
		theme_classic() +
		ylim(-5,5) +
		scale_fill_manual(values=colorRampPalette(c('lightblue','darkblue'))(4)) +
		geom_hline(yintercept=0, linetype='dashed') +
		xlab('') +
		ylab('Relative expression')
	
	print(p)
}	

quartz(w=4.5, h=3)
plotFC(rbind(add_column(rltr.df, subf='RLTR13D5'), add_column(rmer.df, subf='RMER10A')))
plotFC(rbind(add_column(iapc.df, subf='Common'), add_column(iapp.df, subf='Placenta')))
plotFC(rbind(add_column(l1c.df, subf='Common'), add_column(l1p.df, subf='Placenta')))


##ANOVA

summary(aov(plac-tsc ~ dgroup, data=rltr.df))
summary(aov(plac-tsc ~ dgroup, data=rmer.df))
summary(aov(plac-tsc ~ dgroup, data=iapc.df))
summary(aov(plac-tsc ~ dgroup, data=iapp.df))
summary(aov(plac-tsc ~ dgroup, data=l1c.df))
summary(aov(plac-tsc ~ dgroup, data=l1p.df))


##number of genes in each group

group_by(rltr.df, dgroup) %>% summarise(n=length(start))
group_by(rmer.df, dgroup) %>% summarise(n=length(start))
group_by(iapc.df, dgroup) %>% summarise(n=length(start))
group_by(iapp.df, dgroup) %>% summarise(n=length(start))
group_by(l1c.df, dgroup) %>% summarise(n=length(start))
group_by(l1p.df, dgroup) %>% summarise(n=length(start))



############
## ENCODE ##
############

enc = read_tsv('../encode_raw_counts.txt.gz')
enc$Placenta = rowMeans(log2(enc[,grep('placenta',colnames(enc))] + 0.05))
enc$Cerebellum = rowMeans(log2(enc[,grep('cerebellum',colnames(enc))] + 0.05))
enc$Cortex = rowMeans(log2(enc[,grep('cortex',colnames(enc))] + 0.05))
enc$Heart = rowMeans(log2(enc[,grep('heart',colnames(enc))] + 0.05))
enc$Kidney = rowMeans(log2(enc[,grep('kidney',colnames(enc))] + 0.05))
enc$Liver = rowMeans(log2(enc[,grep('liver',colnames(enc))] + 0.05))
enc$SmallIntestine = rowMeans(log2(enc[,grep('smallintestine',colnames(enc))] + 0.05))
enc$Spleen = rowMeans(log2(enc[,grep('spleen',colnames(enc))] + 0.05))
enc$Testis = rowMeans(log2(enc[,grep('testis',colnames(enc))] + 0.05))
enc$Thymus = rowMeans(log2(enc[,grep('thymus',colnames(enc))] + 0.05))

enc2 = select(enc, c(Feature, Placenta:Thymus)) %>%
	pivot_longer(!Feature, names_to='sample', values_to='log2rpm')

##adjust skewed distributions
q75 = group_by(enc2,sample) %>% summarise(quant=quantile(log2rpm,0.75))
enc3 = mutate(enc2, renorm = log2rpm - q75$quant[match(sample, q75$sample)] + mean(q75$quant))


plot.enc = function(data, df, cutoff=1e4) {
	sub = filter(data, Feature %in% df$gene[df$distance<cutoff])
	p = ggplot(sub, aes(x=sample, y=renorm)) +
			geom_boxplot(fill='lightblue', outlier.size=0.5) +
			theme_classic() +
			xlab('') +
			ylab('Log2 RPM') +
			scale_x_discrete(guide = guide_axis(angle = 90))
	
	print(p)
	
	pval = numeric(length(unique(sub$sample))-1)
	i = 1
	for (t in unique(sub$sample)) {
		if (t=='Placenta') next
		test.data = filter(sub, sample=='Placenta' | sample==t)
		pval[i] = wilcox.test(log2rpm~sample, data=test.data, paired=T)$p.value
		i = i+1
	}
	return(p.adjust(pval, method='BH'))
}	

quartz(w=4, h=3)
plot.enc(enc3, rltr.df)
plot.enc(enc3, rmer.df)
plot.enc(enc3, iapp.df)
plot.enc(enc3, iapc.df)
plot.enc(enc3, l1p.df)
plot.enc(enc3, l1c.df)



#############
## Vignola ##
#############

vig = read_tsv('../vignola_raw_counts.txt.gz')
vig$gfp = rowMeans(log2(vig[,grep('_GFP.bam',colnames(vig))] + 0.05))
tm = rowMeans(log2(vig[,grep('_TM.bam',colnames(vig))] + 0.05))
vig$tm = tm + quantile(vig$gfp,0.75) - quantile(tm,0.75)  ##adjust skewed distributions

vig2 = select(vig, c(Feature, ExMesodem=gfp, Trophoblast=tm)) %>%
	pivot_longer(!Feature, names_to='sample', values_to='log2rpm')


plot.vig = function(data, df, cutoff=1e4) {
	sub = filter(data, Feature %in% df$gene[df$distance<cutoff])	
	p =	ggplot(sub, aes(x=sample, y=log2rpm)) +
			geom_boxplot(fill='lightblue', outlier.size=0.5) +
			theme_classic() +
			xlab('') +
			ylab('Log2 RPM') +
			scale_x_discrete(guide = guide_axis(angle = 90))
	
	print(p)
	return(wilcox.test(log2rpm~sample, data=sub, paired=T))
}	

quartz(w=2, h=3)
plot.vig(vig2, rltr.df)
plot.vig(vig2, rmer.df)
plot.vig(vig2, iapp.df)
plot.vig(vig2, iapc.df)
plot.vig(vig2, l1p.df)
plot.vig(vig2, l1c.df)



