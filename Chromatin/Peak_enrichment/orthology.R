library(tidyverse)
library(gplots)


##files from permutation analysis

o.files = list.files(path='orthology',full.names=TRUE)


##read into tibble

over = tibble()
for (f in o.files) {
	ov = read_tsv(f) %>%
		add_column(sample = sub('_overlap.txt','',basename(f))) %>%
		mutate(sig = log.enrich>1 & pval<0.05 & real.overlap>10,
			mark = strsplit(sample, '_')[[1]][1],
			type = strsplit(sample, '_')[[1]][2],
			species = strsplit(sample, '_')[[1]][3],)
	over = rbind(over, ov)
}
colnames(over)[1] = 'rep'


##filter significantly enriched families from orthologous groups

epi = filter(over, type=='epiOrthologous')
gen = filter(over, type=='genOrthologous')

epi.sig = filter(epi, rep %in% unique(epi$rep[epi$sig]))
gen.sig = filter(gen, rep %in% unique(gen$rep[gen$sig])) %>%	
		filter(rep!='RLTR26B_MM') #unnanotated in pahari


##heatmaps

enr.hm = function(data) {
	enr.sum = select(data, rep, log.enrich, sample) %>%
		pivot_wider(names_from=sample, values_from=log.enrich)
	
	enr.mat = as.matrix(enr.sum[,-1])
	row.names(enr.mat) = enr.sum$rep
	enr.mat[enr.mat==-Inf] = min(enr.mat[enr.mat!=-Inf])
	
	heatmap.2(enr.mat, scale='none',
			dendrogram='row', Colv=FALSE,
			col=colorRampPalette(c('navy','black','yellow'))(100),
			breaks=seq(-4,4,length.out=101),
			trace='none',
			margins=c(10,7),
			density.info='none',
			key.title='log2 enrichment',
			cexCol=0.8)
}


quartz(h=3,w=3)
enr.hm(epi.sig)

quartz(h=4,w=3)
enr.hm(gen.sig)


##non-orthologous elements

mm = filter(over, type=='nonOrthologous', species=='mm')
#mp = filter(over, type=='nonOrthologous', species=='mp')

sel.te = c('IAPLTR1_Mm','IAPLTR2_Mm','IAPLTR2a','L1MdTf_I','L1MdTf_II', 'L1MdA_I')

filter(mm, rep %in% sel.te) %>%
ggplot(aes(x=rep, y=real.overlap/random.mean, fill=mark)) +
	geom_bar(stat='identity', position=position_dodge(), col='black') +
	theme_classic() +
	xlab('') +
	ylab('observed/expected') +
	geom_hline(yintercept=1, linetype='dashed')




##RMER10A barplot

rmer = filter(over, rep=='RMER10A', mark=='H3K27ac', type!='nonOrthologous') %>%
	mutate(enrich=real.overlap/random.mean)

ggplot(rmer, aes(x=type, y=enrich, fill=species)) +
	geom_bar(stat='identity', position=position_dodge(), col='black') +
	theme_classic() +
	xlab('') +
	ylab('observed/expected') +
	geom_hline(yintercept=1, linetype='dashed')


