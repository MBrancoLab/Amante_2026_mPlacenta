library(tidyverse)
library(gplots)


##read IAP data
##rpm counts from Seqmonk using IAP_merged annotation, only intergenic elements

iap = read_tsv('IAP_intergenic_RPM.txt')
iap$placenta = rowMeans(select(iap, starts_with('Placenta')))


##mark loci with minimum expression and DE

iap = mutate(iap, expr = placenta>-3 | GFPTSCs_musculus.bam>-3 | bTSCs_musculus.bam>-3,
			up = placenta-GFPTSCs_musculus.bam>1 & placenta-bTSCs_musculus.bam>1,
			down = GFPTSCs_musculus.bam-placenta>1 & bTSCs_musculus.bam-placenta>1)


##plot heatmap

iap.mat = as.matrix(filter(iap, expr) %>% select(ends_with('bam')))

heatmap.2(iap.mat, trace='none', scale='none',
	col=hcl.colors(100), density.info='none')


##count subfamilies and conformation

subf.count = group_by(iap, Feature) %>%
	summarise(exp=sum(expr), up=sum(expr&up), down=sum(expr&down))

conf.count = group_by(iap, Description) %>%
	summarise(exp=sum(expr), up=sum(expr&up), down=sum(expr&down))


##plot subfamily counts

filter(subf.count, !grepl('int',Feature) & !grepl('_I',Feature)) %>%
pivot_longer(cols=!Feature, names_to='type', values_to='counts') %>%
ggplot(aes(x=Feature, y=counts, fill=factor(type, levels=c('exp','up','down')))) +
	geom_bar(stat='identity', position=position_dodge(), col='black') +
	theme_classic() +
	xlab('') +
	ylab('Number of elements') +
	labs(fill='') +
	scale_x_discrete(guide = guide_axis(angle = 90)) +
	scale_fill_manual(values=c('grey','orange','lightblue'))


##plot element lengths

mutate(iap, len=End-Start) %>% filter(expr, len<10000) %>%
ggplot(aes(x=len, fill=Description)) +
	geom_histogram(col='black') +
	theme_classic() +
	xlab('Length (bp)') +
	ylab('Number of elements') +
	scale_fill_manual(values=c('black','orange','lightblue','grey'))

