library(tidyverse)


##stringtie quantification data

q.files = list.files('merged_quant', full.names=TRUE)

quant = tibble()
for (f in q.files) {
	sub = read_tsv(f) %>%
		add_column(sample = gsub('_quant.ctab.gz','',basename(f)))%>%
		mutate(tissue = gsub('-[[:digit:]]','',sample))
	quant = rbind(quant,sub)
}



##select for genes with IAP chimeras
##additional info added manually afterwards

#sel.tx = read_tsv('iap_l1_transcripts_curated.txt') %>% filter(grepl('IAP', te.name))
#gene.ids = unique(quant$gene_id[quant$gene_name %in% sel.tx$gene.name])
#sub.quant = filter(quant, tissue=='bTSCs_musculus', gene_id %in% gene.ids) %>%
#		select(t_id:gene_name)
#write_tsv(sub.quant, 'iap_chimera_isoforms.txt')


##add info on which isoforms from selected genes have IAP TSS
##calculate percentage of IAP TSS usage

meta = read_tsv('iap_chimera_isoforms.txt')

iap = inner_join(select(quant, !gene_name), select(meta, t_id, gene_name, te_tss), by='t_id')

iap.sum = group_by(iap, gene_name, sample, tissue) %>%
	summarise(fpkm = log2(sum(FPKM)+0.1),
			perc = case_when(sum(FPKM)==0 ~ 0,
				.default = sum(FPKM[te_tss])/sum(FPKM)*100)) %>%	
	group_by(gene_name, tissue) %>%
	summarise(fpkm=mean(fpkm), perc=mean(perc)) %>%
	pivot_wider(names_from=tissue, values_from=c(fpkm,perc))

ggplot(iap.sum, aes(x=perc_Placenta_musculus-perc_bTSCs_musculus, y=fpkm_Placenta_musculus-fpkm_bTSCs_musculus)) +
	geom_point(col='navy') +
	geom_hline(yintercept=0, linetype='dashed') +
	theme_classic() +
	xlab('Difference in IAP usage (%)') +
	ylab('Placenta/TSC log2 ratio')



