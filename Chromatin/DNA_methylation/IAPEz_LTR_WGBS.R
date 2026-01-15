library(tidyverse)

iap = read_tsv('IAPEz_LTR_WGBS.txt')


ggplot(arrange(iap,`Feature Orientation`), aes(x=Brain, y=Placenta)) +
	geom_point(aes(col=`Feature Orientation`, size=`Feature Orientation`)) +
	geom_abline(slope=1, intercept=0, linetype='dashed') +
	theme_classic() +
	scale_color_manual(values=c('grey','orange'), name='H3K4me3') +
	scale_size_manual(values=c(0.3, 1), name=NULL)

ggplot(iap, aes(x='Placenta', y=Placenta)) +
	geom_violin() +
	geom_point(aes(x='Placenta', y=Placenta[Type=='H3K4me3']), position=position_jitter(0.2))

