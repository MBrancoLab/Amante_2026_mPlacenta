library(pzfx)
library(tidyverse)


##read data from prism file

plac = 'VM-IAP_tail_placenta_methylation.pzfx'
ptables = pzfx_tables(plac)
col_tables = ptables[grep('column',ptables)]

met = tibble()
for (tab in col_tables) {
	df = read_pzfx(plac, tab)
	colnames(df)[1] = 'ID'
	df.long = pivot_longer(df, cols=!ID, names_to='tissue', values_to='methylation') %>%
		add_column(gene = sub('_avg_column','',tab)) %>%
		mutate(tlabel = case_when(grepl('Placenta',tissue) ~ 'P',
								grepl('Tail',tissue) ~ 'T'))
	met = rbind(met, df.long)
}


##plot

met_plot = function(data) {
	ggplot(data, aes(x=factor(tlabel, levels=c('T','P')), y=methylation, colour=factor(tlabel, levels=c('T','P')))) +
		geom_line(aes(group = interaction(ID, gene)), colour = "grey", linetype='dashed') +
		geom_point(position = position_dodge(0.7)) +
		ylim(0, 100) +
		facet_wrap(vars(gene), nrow=1) +
		xlab('') +
		ylab('% Methylation') +
		theme_classic() +
		scale_colour_manual(values = c('black', 'orange'), name = NULL)
}

met_plot(filter(met, gene %in% c('Tfpi','Mbln1','Bmf','Marveld2','Gm13849','Rab6b','Fam78b','Eps8l1','Rnf157')))
met_plot(filter(met, gene %in% c('Cdc73','H2T24','Diap3_nr','Dst','Cdk15','Ell2','Slc24a3','Zak')))


##paired t-tests

p = numeric()

for (i in 1:length(col_tables)) {
	df = read_pzfx(plac, col_tables[i])
	colnames(df)[1] = 'ID'
	p[i] = t.test(df[,2],df[,3], paired=TRUE)$p.value
}
names(p) = sub('_avg_column','',col_tables)

padj = p.adjust(p, method='BH')