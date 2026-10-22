library(tidyverse)


##Suv39h2 isoforms

suv = read_tsv('Suv39h2_qPCR.txt')

ggplot(suv, aes(x=Isoform, y=Expression, fill=Sample)) +
	geom_bar(position=position_dodge(), stat='summary', fun=mean, col='black', width=0.6) +
	geom_point(position=position_dodge(width=0.6), col='black') +
	xlab('') +
	ylab('Relative expression') +
	theme_classic()


##TSC differentiation

diff = read_tsv('diffTSC_qPCR.txt')

quartz(w=3.5, h=2)
filter(diff, Target %in% c('Cdx2','Elf5','Eomes')) %>%
ggplot(aes(x=Sample, y=Target, size=Normalised)) +
	geom_point() +
	scale_size_continuous(range = c(0,7)) +
	xlab('') +
	ylab('') +
	theme_classic()

quartz(w=3.5, h=2.5)
filter(diff, Target %in% c('Prl3d1','Syna','Synb','Ctsq')) %>%
ggplot(aes(x=Sample, y=Target, size=Normalised)) +
	geom_point() +
	scale_size_continuous(range = c(0,7)) +
	xlab('') +
	ylab('') +
	theme_classic()