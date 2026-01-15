library(tidyverse)


##sgRNA hits predicted with bowtie2 (-a option)
##bam converted to bed and intersected with repeatmasker (-wao option)

hits = read_tsv('sgRNA_rmsk_hits.bed', col_names=F)


##get target frequencies

freq = group_by(hits, X10) %>%
	summarise(count = length(X10)) %>% 
	mutate(freq = count/nrow(hits))


##pool low frequency targets

lowf = filter(freq, count<50)

freq2 = rbind(filter(freq, count>50),
	tibble(X10='Others', count=sum(lowf$count), freq=sum(lowf$freq)))


##plot

ggplot(freq2, aes(x=2, y=freq, fill=X10)) +
	geom_bar(stat='identity') +
	coord_polar(theta = "y", start = 0) +
	xlim(0.5, 2.5) +
	theme_void()