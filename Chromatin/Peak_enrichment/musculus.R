library(tidyverse)
library(gplots)


##files from permutation analysis

o.files = c('musculus_cuttag/H3K27ac_bTSCs_overlap.txt',
			'musculus_cuttag/H3K27ac_GFPTSCs_overlap.txt',
			'musculus_cuttag/H3K27ac_placenta_overlap.txt',
			'musculus_cuttag/H3K4me1_bTSCs_overlap.txt',
			'musculus_cuttag/H3K4me1_GFPTSCs_overlap.txt',
			'musculus_cuttag/H3K4me1_placenta_overlap.txt',
			'musculus_cuttag/H3K4me1_trophoblast_overlap.txt',
			'musculus_cuttag/H3K4me1_endothelial_overlap.txt',
			'musculus_cuttag/H3K4me3_bTSCs_overlap.txt',
			'musculus_cuttag/H3K4me3_GFPTSCs_overlap.txt',
			'musculus_cuttag/H3K4me3_placenta_overlap.txt',
			'musculus_cuttag/H3K4me3_trophoblast_overlap.txt',
			'musculus_cuttag/H3K4me3_endothelial_overlap.txt',
			'musculus_atac/ATAC_TS_overlap.txt',
			'musculus_atac/ATAC_d2diff_overlap.txt',
			'musculus_atac/ATAC_trophoblast_overlap.txt',
			'musculus_atac/ATAC_endothelial_overlap.txt')


##read into tibble

over = tibble()
for (f in o.files) {
	ov = read_tsv(f) %>%
		add_column(sample = gsub('_overlap.txt','',basename(f))) %>%
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

quartz(h=5.7,w=4.6)
enr.hm(filter(sig, grepl('H3K27ac',sample) | grepl('ATAC',sample)))

quartz(h=5.7,w=4.2)
enr.hm(filter(sig, grepl('H3K4me1',sample)))
enr.hm(filter(sig, grepl('H3K4me3',sample)))


##enrichment scatter plots

enr.scatter = function(sample_to_plot, ymax=35) {
	filter(over, sample==sample_to_plot) %>% arrange(sig) %>%
	ggplot(aes(x=log10(real.overlap), y=real.overlap/random.mean)) +
		geom_point(aes(colour=sig, size=sig)) +
		theme_classic() +
		scale_color_manual(values=c('grey','black')) +
		scale_size_manual(values=c(0.5,1)) +
		ylim(0,ymax) +
		xlab('log10 n overlapped') +
		ylab('observed/expected')
}

quartz(w=4, h=2.7)
enr.scatter('H3K27ac_GFPTSCs')
enr.scatter('H3K27ac_placenta')

enr.scatter('H3K27ac_bTSCs')
enr.scatter('H3K4me1_placenta')
enr.scatter('H3K4me3_placenta')


