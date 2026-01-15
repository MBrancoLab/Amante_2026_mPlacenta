##Plot average normalised read count across selected TE families
##Source data produced by deeptools (trends.sh script)

library(tidyverse)
library(gplots)


##read data

trend_files = list.files('Trends', full.names=TRUE)

data = tibble()
for (f in trend_files) {
	sub = read_tsv(f, skip=2, col_names=c('sample','ignore',1:100)) %>%
		select(!ignore) %>%
		pivot_longer(cols=!sample,names_to='bin') %>%
		mutate(across(bin, as.integer)) %>%
		add_column(te = gsub('_trend.txt','',basename(f)))
	data = rbind(data, sub)
}


##plot trends

plot_trend = function(data, subf, samples, ymax=max(data$value)) {
	filter(data, te==subf, sample %in% samples) %>%
	ggplot(aes(x=bin, y=value, colour=sample)) +
		geom_line(linewidth=1) +
		ylim(0, ymax) +
		theme_classic() +
		ylab('Mean RPM') +
		scale_x_continuous(name='', breaks=c(0,20,80,100), labels=c('-2kb','start','end','+2kb'))
}


##plots

k4_multi = c('H3K4me3_bTSCs_musculus_multi','H3K4me3_GFPTSCs_musculus_multi',
			'H3K4me3_placenta_musculus_multi','H3K4me3_trophoblast_musculus_multi',
			'H3K4me3_endothelial_musculus_multi')

k4_unique = gsub('_multi','',k4_multi)

plot_trend(data, subf='RLTR4', samples=k4_multi, ymax=4)
plot_trend(data, subf='RLTR4', samples=k4_unique, ymax=4)

plot_trend(data, subf='L1MdTf_I', samples=k4_multi, ymax=1)
plot_trend(data, subf='L1MdTf_I', samples=k4_unique, ymax=1)

plot_trend(data, subf='IAPLTR1_Mm', samples=k4_multi, ymax=1)
plot_trend(data, subf='IAPLTR1_Mm', samples=k4_unique, ymax=1)


##extract 5' UTR enrichment

utr = filter(data, bin>=20 & bin<=80) %>%
	group_by(sample, te) %>% summarise(summit = max(value))
	

##plot heatmap

sel = c('H3K4me3_bTSCs_musculus_multi',
		'H3K4me3_GFPTSCs_musculus_multi',
		'H3K4me3_placenta_musculus_multi',
		'H3K4me3_trophoblast_musculus_multi',
		'H3K4me3_endothelial_musculus_multi')

utr.sel = filter(utr, sample %in% sel) %>%
	pivot_wider(names_from=te, values_from=summit)

utr.mat = as.matrix(utr.sel[,-1])
rownames(utr.mat) = utr.sel$sample

heatmap.2(utr.mat[sel,], scale='none',
		Rowv=FALSE, dendrogram='none',
		col = hcl.colors(60),
		breaks = seq(0,2,length.out=61),
		margins = c(8,8), cexRow=0.7, cexCol=0.7,
		trace = 'none', density.info='none')

