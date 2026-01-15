library(tidyverse)


##read squire files

sq = tibble()
for (f in list.files(path='./squire')) {
	sub = read_tsv(paste('./squire',f,sep='/')) %>%
		add_column(group=sub('_0[[:digit:]]_subFcounts.txt','',f))
	sq = rbind(sq, sub)
}
colnames(sq)[3] = 'te'


##plot selected IAP families

sel.te = c('IAPLTR1_Mm:ERVK:LTR', 'IAPLTR1a_Mm:ERVK:LTR',
	'IAPLTR2_Mm:ERVK:LTR', 'IAPLTR2a:ERVK:LTR',
	'IAPLTR2a2_Mm:ERVK:LTR', 'IAPLTR2b:ERVK:LTR')

sq.sel = filter(sq, te %in% sel.te)

ggplot(sq.sel, aes(x=te, y=fpkm, fill=factor(group,levels=c('KRAB','IAP_gRNA')))) +
	geom_bar(stat='summary', fun=mean, col='black', position=position_dodge()) +
	geom_point(position=position_dodge(0.9)) +
	theme_classic() +
	xlab('') +
	ylab('FPKM') +
	labs(fill = '')


##anova

model = with(sq.sel, aov(fpkm ~ group * te))
summary(model)

TukeyHSD(model)
