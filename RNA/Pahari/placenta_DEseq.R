library(tidyverse)
library(DESeq2)
library(topGO)
library(org.Mm.eg.db)


##gene raw counts from Seqmonk

mm = read_tsv('../musculus_raw_counts.txt.gz') %>%
	dplyr::select(!c(bTSCs_musculus, GFPTSCs_musculus))

mp = read_tsv('../pahari_raw_counts.txt.gz') %>%
	dplyr::select(Probe, Placenta_pahari_2, Placenta_pahari_3)
colnames(mp)[1] = 'pahari_id'


##orthology info from biomart

ortho = read_tsv('mm.mp.orthologues.txt.gz')
ortho1 = filter(ortho, `Shrew mouse homology type`=='ortholog_one2one') %>%
	dplyr::select(`Gene stable ID`, `Shrew mouse gene stable ID`) %>%
	distinct()
colnames(ortho1) = c('ID', 'pahari_id')


##merge datasets
##remove sex chromosome data

mm.ortho = inner_join(mm, ortho1, by='ID') %>%
	relocate(pahari_id, .after=ID)

mm.mp = inner_join(mm.ortho, mp, by='pahari_id') %>%
	filter(Chromosome!='X', Chromosome!='Y')


##extract count data for DESeq2

counts = data.frame(dplyr::select(mm.mp, Placenta_musculus_1:Placenta_pahari_3))
rownames(counts) = mm.mp$Feature

coldata = DataFrame(sample=factor(sub('_[[:digit:]]','',colnames(counts))))


##get DE genes

dds = DESeqDataSetFromMatrix(counts, coldata, design=~sample)
dds = DESeq(dds)

res = results(dds, contrast=c('sample','Placenta_musculus','Placenta_pahari'))
de = subset(res, padj<0.05 & abs(log2FoldChange)>1) %>%
	as_tibble(rownames='gene') %>%
	mutate(species = case_when(log2FoldChange>0 ~ 'musculus', log2FoldChange<0 ~ 'pahari')) %>%
	arrange(-log2FoldChange)

group_by(de, species) %>% summarise(n.de=length(gene))

write_tsv(de, 'placenta_DE_genes.txt')


##change in Marsh markers

markers = c('Ror2','Lgr5','Met','Gcm1','Gcgr','Vegfa','Egfr','Tcf7l1','Epha4',
	'Tgfa','Glis1','Stra6','Nos1ap','Podxl','Lepr','Cdh4','Prune2',
	'Ncam1','Igfbp7','Pla2g4d','Plac8','Mitf','Flt1','Prl8a9','Slco2a1')

marsh = as_tibble(res[rownames(res) %in% markers,], rownames='gene')

ggplot(marsh, aes(y=factor(gene, levels=markers), x=log2FoldChange)) +
	geom_bar(stat='identity', colour='black', fill='grey') +
	geom_vline(xintercept=0) +
	ylab('') +
	theme_classic()

cpm = as_tibble(fpm(dds), rownames='gene') %>%
	pivot_longer(cols=!gene, names_to='sample', values_to='expr') %>%
	mutate(species = sapply(strsplit(sample,split='_'), function(x) x[2]))

plot_markers = function(markers) {
	filter(cpm, gene %in% markers) %>%
	ggplot(aes(y=factor(gene, levels=markers), x=log2(expr), fill=species)) +
		geom_bar(stat='summary', fun=mean, colour='black', position=position_dodge()) +
		geom_point(position=position_dodge(0.9), size=0.5) +
		theme_classic() +
		ylab('') +
		xlab('log2 CPM')
}

plot_markers(c('Slco2a1','Prl8a9','Pla2g4d','Cdh4'))
plot_markers(c('Gcgr','Vegfa','Met','Ror2'))


##gene ontology analysis

go.ann = function(de.genes) {
	genes = factor(as.integer(rownames(counts) %in% de.genes$gene))
	names(genes) = rownames(counts)
	go = new('topGOdata', ontology = 'BP',
		allGenes = genes, geneSelectionFun = function(x) x==0,
		annotationFun = annFUN.org, mapping = 'org.Mm.eg.db', ID = 'SYMBOL')
	return(go)
}

go.plot = function(go, nodes=5) {
	rt = runTest(go, algorithm='classic', statistic='fisher')
	gt = GenTable(go, classicFisher=rt, topNodes=nodes)

	ggplot(gt, aes(x=-log10(as.numeric(classicFisher)), y=Term)) +
		geom_bar(stat='identity', fill='grey', colour='black', width=0.7) +
		theme_classic() +
		ylab('') +
		xlab('-log10 p')
}

go.mm = go.ann(filter(de, species=='musculus'))
go.plot(go.mm)

go.mp = go.ann(filter(de, species=='pahari'))
go.plot(go.mp)

