library(tidyverse)
library(gplots)


##files from permutation analysis

o.files = list.files(path='encode', full.names=TRUE)


##read into tibble

over = tibble()
for (f in o.files) {
	ov = read_tsv(f) %>%
		add_column(sample = gsub('_overlap.txt','',basename(f))) %>%
		mutate(sig = log.enrich>1 & pval<0.05 & real.overlap>10,
			tissue = sub('H3K27ac_', '', sample))
	over = rbind(over, ov)
}
colnames(over)[1] = 'rep'


##plot significantly enriched subfamilies in each tissue

t.levels = c('bonemarrow','heart','kidney','liver','smallintestine','spleen','thymus',
		'olfactorybulb','cortex','cerebellum','testis','placenta','hplacenta')

ggplot(filter(over, sig), aes(x=factor(tissue, levels=t.levels), y=log.enrich, col=sample)) +
	geom_point(position=position_jitter(0.2), show.legend=FALSE) +
	theme_classic() +
	xlab('') +
	ylab('log2 observed/expected') +
	scale_x_discrete(guide = guide_axis(angle = 90))


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
enr.scatter('H3K27ac_placenta')
enr.scatter('H3K27ac_testis')
enr.scatter('H3K27ac_hplacenta')


##compare enrichments between ENCODE and own data

comp = full_join(read_tsv('musculus_cuttag/H3K27ac_placenta_overlap.txt'),
		read_tsv('encode/H3K27ac_placenta_overlap.txt'), by='...1')
colnames(comp)[1] = 'rep'

filter(comp, real.overlap.x>10, real.overlap.y>10) %>%
ggplot(aes(x=log.enrich.x, y=log.enrich.y)) +
	geom_point(col='grey', size=1) +
	theme_classic() +
	xlab('This study') +
	ylab('ENCODE')



