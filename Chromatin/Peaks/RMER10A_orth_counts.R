library(tidyverse)

##Read counts

mm = read_tsv('RMER10A_orthMM_counts.txt') %>%
	mutate(Peak_mm = `Feature Orientation`=='overlapping') %>%
	select(Probe:End,Peak_mm,H3K27ac_placenta_musculus.bam)
	
mp = read_tsv('RMER10A_orthMP_counts.txt') %>%
	mutate(Peak_mp = `Feature Orientation`=='overlapping') %>%
	select(Probe:End,Peak_mp,H3K27ac_placenta)

counts = inner_join(mm, mp, by='Probe') %>%
	mutate(any_peak = Peak_mm | Peak_mp)


##Plot

filter(counts, any_peak) %>%
ggplot(aes(x=H3K27ac_placenta_musculus.bam, y=H3K27ac_placenta)) +
	geom_point(colour='navy') +
	theme_classic() +
	xlab('H3K27ac musculus (RPM)') +
	ylab('H3K27ac pahari (RPM)')


