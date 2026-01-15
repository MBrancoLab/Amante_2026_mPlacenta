library(tidyverse)


##read IAP expression data
##rpm counts from Seqmonk using IAP_merged annotation, only intergenic elements

iap = read_tsv('IAP_intergenic_RPM.txt')
iap$placenta = rowMeans(select(iap, starts_with('Placenta')))


##select well expressed proviral loci

sel = filter(iap, placenta>1 | GFPTSCs_musculus.bam>1 | bTSCs_musculus.bam>1,
			Description=='proviral',
			grepl('IAPLTR1',Feature) | grepl('IAPLTR2',Feature) | grepl('IAPEY4',Feature)) %>%
	mutate(id = paste(Feature,'_chr',Chromosome,'_',Start,sep=''))


##parse into gtf format

gtf = data.frame(chr = paste('chr',sel$Chromosome,sep=''),
	source = rep('unknown',nrow(sel)),
	feature = rep('exon',nrow(sel)),
	start = sel$Start,
	end = sel$End,
	score = rep('.',nrow(sel)),
	strand = sel$`Feature Strand`,
	frame = rep('.',nrow(sel)),
	attribute = paste('gene_id "',sel$id,'"; transcript_id "',sel$id,'"',sep=''))


##write

write_tsv(gtf, 'expressed_IAPs.gtf', col_names=FALSE, escape='none')