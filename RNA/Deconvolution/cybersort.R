
library(tidyverse)


##Cell type fractions from CYBERSORT output
##TO ADD: some details on cybersort processing

cs = read_tsv('cybersort.txt')


##establish tissue order and colours

torder = c('LaTP', 'SynTI', 'SynTII', 'S-TGC', 'JZP', 'SpT', 'glycogen cell',
			'pro-angiogenic endothelial', 'lymphogenic endothelial', 'decidual stroma', 'fetal mesenchyme',
			'B-cell', 'T-cell', 'natural killer cell', 'macrophage', 'erythrocyte')

tcol = c('#282e30', '#7e7f7b', '#b2b5af', '#ebebe9', '#c50000', '#ec2d2e', '#ff6347',
		'#1e3f7a', '#6894c7', '#006400', '#7fe240',
		'#ba741d', '#c75310', '#fda840', '#ecc729', '#f9eb6c')


##get means

cs.av = group_by(cs, tissue, celltype) %>% summarise(av.frac = mean(fraction))


##plot

filter(cs.av, tissue!='TSCs') %>%
ggplot(aes(x=tissue, y=av.frac, fill=factor(celltype, levels=torder))) +
	geom_bar(stat='identity') +
	scale_fill_manual(values=tcol) +
	xlab('') +
	ylab('') +
	theme_classic()

