library(tidyverse)
library(patchwork)


##read squire files

sq.mus = tibble()
for (f in list.files('musculus', pattern='subFcounts.txt')) {
	sub = read_tsv(paste('musculus', f, sep='/')) %>%
		mutate(tissue=sub('_musculus[[:print:]]+','',Sample))
	sq.mus = rbind(sq.mus, sub)
}

sq.enc = tibble()
for (f in list.files('encode', pattern='subFcounts.txt')) {
	sub = read_tsv(paste('encode', f, sep='/')) %>%
		mutate(tissue=sub('-[[:digit:]]','',sub('RNA_','',Sample)))
	sq.enc = rbind(sq.enc, sub)
}

sq.vig = tibble()
for (f in list.files('vignola', pattern='subFcounts.txt')) {
	sub = read_tsv(paste('vignola', f, sep='/')) %>%
		mutate(tissue = case_when(grepl('GFP',Sample) ~ 'endothelial',
								grepl('TM',Sample) ~ 'trophoblast'))
	sq.vig = rbind(sq.vig, sub)
}


##read TEtranscripts files

te.mus = tibble()
for (f in list.files('musculus', pattern='TEtranscripts.cntTable')) {
	sub = read_tsv(paste('musculus', f, sep='/'), col_names=c('roi','count'), skip=1) %>%
		mutate(rpm=count/sum(count)*1e6, tissue=sub('_musculus[[:print:]]+','',f)) %>%
		add_column(sample=gsub('_TEtranscripts.cntTable','',f))
	te.mus = rbind(te.mus, sub)
}

te.enc = tibble()
for (f in list.files('encode', pattern='TEtranscripts.cntTable')) {
	sub = read_tsv(paste('encode', f, sep='/'), col_names=c('roi','count'), skip=1) %>%
		mutate(rpm=count/sum(count)*1e6, tissue=sub('-[[:print:]]+','',sub('RNA_','',f))) %>%
		add_column(sample=gsub('_TEtranscripts.cntTable','',f))
	te.enc = rbind(te.enc, sub)
}

te.vig = tibble()
for (f in list.files('vignola', pattern='TEtranscripts.cntTable')) {
	sub = read_tsv(paste('vignola', f, sep='/'), col_names=c('roi','count'), skip=1) %>%
		mutate(rpm=count/sum(count)*1e6,
			tissue=case_when(grepl('GFP',f) ~ 'endothelial',
						grepl('TM',f) ~ 'trophoblast')) %>%
		add_column(sample=gsub('_TEtranscripts.cntTable','',f))
	te.vig = rbind(te.vig, sub)
}


##plot individual TEs

plot.te = function(te, squire, tetr) {
	p1 = squire %>% filter(`Subfamily:Family:Class`==te) %>%
	ggplot(aes(x=tissue, y=fpkm)) +
		geom_bar(stat='summary', fun=mean) +
		geom_point() +
		theme_classic() +
		xlab('') +
		ylab('FPKM') +
		ggtitle('SQuIRE')
	
	p2 = tetr %>% filter(roi==te) %>%
	ggplot(aes(x=tissue, y=rpm)) +
		geom_bar(stat='summary', fun=mean) +
		geom_point() +
		theme_classic() +
		xlab('') +
		ylab('RPM') +
		ggtitle('TEtranscripts')
	
	p1 + p2
}


plot.te('IAPLTR1_Mm:ERVK:LTR', sq.mus, te.mus)
plot.te('IAPLTR1_Mm:ERVK:LTR', sq.enc, te.enc)
plot.te('IAPLTR1_Mm:ERVK:LTR', sq.vig, te.vig)

plot.te('L1Md_T:L1:LINE', sq.mus, te.mus)
plot.te('L1Md_T:L1:LINE', sq.enc, te.enc)
plot.te('L1Md_T:L1:LINE', sq.vig, te.vig)

