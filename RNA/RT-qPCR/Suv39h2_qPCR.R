library(tidyverse)


##read data

suv = read_tsv('Suv39h2_qPCR.txt')


##plot

ggplot(suv, aes(x=Isoform, y=Expression, fill=Sample)) +
	geom_bar(position=position_dodge(), stat='summary', fun=mean, col='black', width=0.6) +
	geom_point(position=position_dodge(width=0.6), col='black') +
	xlab('') +
	ylab('Relative expression') +
	theme_classic()
