library(tidyverse)
library(gplots)


##files from permutation analysis

o.files = list.files(path='pahari_cuttag',full.names=TRUE)


##read into tibble

over = tibble()
for (f in o.files) {
	ov = read_tsv(f) %>%
		add_column(sample = gsub('_pahari_overlap.txt','',basename(f))) %>%
		mutate(sig = log.enrich>1 & pval<0.05 & real.overlap>10)
	over = rbind(over, ov)
}
colnames(over)[1] = 'rep'


##filter significantly enriched families in at least one sample

sig = filter(over, rep %in% unique(over$rep[over$sig]))


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
			key.title='log2 enrichment')
}

quartz(h=4.5,w=4.5)
enr.hm(sig)


