library(tidyverse)


##counts from liftover.sh script

counts = read_tsv('orthology_counts.txt')


##plot

sp.order = c('rat','pahari','caroli','musculus')
te.order = c('RMER10A','RLTR13D5','IAP','L1MdA','L1MdTf','L1MdGf')

ggplot(counts, aes(x=factor(te,levels=te.order), y=factor(species,levels=sp.order))) +
	geom_point(aes(size=percentage),stroke=0) +
	scale_size_continuous(range=c(0,8)) +
	scale_x_discrete(guide=guide_axis(angle=90)) +
	theme_classic() +
	xlab('') +
	ylab('')