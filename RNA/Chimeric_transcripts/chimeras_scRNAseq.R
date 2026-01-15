library(Seurat)
library(tidyverse)



##Load Seurat object

plac = LoadSeuratRds('../scRNA_processing/Marsh_IAP_Seurat.Rds')


##Chimeras

chim = read_tsv('iap_l1_transcripts_curated.txt')
cgenes = unique(chim$gene.name)


##chimera expression summary

plot.means = function(data, genes) {
	idents = Idents(data)
	sub = as_tibble(as.matrix(plac[['RNA']]$counts), rownames='gene') %>%
		filter(gene %in% genes) %>%
		pivot_longer(cols=!gene, names_to='barcode', values_to='expr') %>%
		mutate(cluster = idents[match(barcode, names(idents))])

	expr.sum = group_by(sub, gene, cluster) %>%
		summarise(av_expr=mean(expr), perc_cells=sum(expr>1)/length(expr)*100)

	c.order = c("SpT", "Glycogen cells", "S-TGC", "SynTI", "SynTII",
	"Endothelial", "Pericytes", "Mesenchymal I", "Mesenchymal II", "Decidua")
	
	g.order = filter(expr.sum, cluster=='SpT') %>% arrange(perc_cells) %>% pull(gene)

	ggplot(expr.sum, aes(x=factor(gene,levels=g.order), y=factor(cluster,levels=c.order))) +
		geom_point(aes(size=perc_cells, colour=av_expr), stroke=0) +
		scale_size_continuous(range = c(0,7)) +
		scale_colour_gradient(low='grey', high='red') +
		theme_classic() +
		theme(axis.text.x = element_text(angle = 45, hjust=1)) +
		xlab('') +
		ylab('')
}

plot.means(plac, 'Noct')

plot.means(plac, cgenes[cgenes!='Noct'])

no.expr = c('4930481A15Rik', 'Csta2', 'Gm20110', 'Gm36888')
plot.means(plac, cgenes[!(cgenes %in% no.expr) & cgenes!='Noct'])





