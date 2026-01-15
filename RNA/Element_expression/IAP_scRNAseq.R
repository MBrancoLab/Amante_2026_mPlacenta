library(Seurat)
library(dplyr)
library(tidyr)
library(gplots)


##Load Seurat object

plac = LoadSeuratRds('../scRNA_processing/Marsh_IAP_Seurat.Rds')


##IAP expression

idents = Idents(plac)
iap = as_tibble(as.matrix(plac[['RNA']]$counts), rownames='gene') %>%
	filter(grepl('IAP', gene)) %>%
	pivot_longer(cols=!gene, names_to='barcode', values_to='expr') %>%
	mutate(cluster = idents[match(barcode, names(idents))])

iap.sum = group_by(iap, gene, cluster) %>%
	summarise(perc_cells=sum(expr>1)/length(expr)*100) %>%
	pivot_wider(names_from=cluster, values_from=perc_cells)

iap.mat = as.matrix(iap.sum[,-1])
rownames(iap.mat) = iap.sum$gene

c.order = c("SpT", "Glycogen cells", "S-TGC", "SynTI", "SynTII",
	"Endothelial", "Pericytes", "Mesenchymal I", "Mesenchymal II", "Decidua")

heatmap.2(iap.mat[order(rowMeans(iap.mat), decreasing=TRUE),c.order],
	scale='none',
	Rowv=FALSE, Colv=FALSE,
	dendrogram='none',
	col=hcl.colors(100),
	trace='none',
	density.info='none',
	margins=c(11,9))

VlnPlot(plac, features = c("IAPLTR1-Mm-chr12-108079292","IAPLTR1-Mm-chr10-78020226","IAPLTR1-Mm-chr5-124104495"), pt.size=0)


