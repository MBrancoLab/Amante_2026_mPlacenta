library(tidyverse)


##normalise read counts (from Seqmonk)

enc = read_tsv('../encode_raw_counts.txt.gz') %>%
	mutate(across(cerebellum_encode_1:thymus_encode_2, ~log2(./sum(.)*1e6+0.05)))

mp = read_tsv('../pahari_raw_counts.txt.gz') %>%
	dplyr::select(Probe, Brain_pahari:Spleen_pahari) %>%
	mutate(across(Brain_pahari:Spleen_pahari, ~log2(./sum(.)*1e6+0.05)))
colnames(mp)[1] = 'pahari_id'


##orthology info from biomart

ortho = read_tsv('mm.mp.orthologues.txt.gz')
ortho1 = filter(ortho, `Shrew mouse homology type`=='ortholog_one2one') %>%
	dplyr::select(`Gene stable ID`, `Shrew mouse gene stable ID`) %>%
	distinct()
colnames(ortho1) = c('ID', 'pahari_id')


##merge datasets
##remove sex chromosome data

enc.ortho = inner_join(enc, ortho1, by='ID') %>%
	relocate(pahari_id, .after=ID)

enc.mp = inner_join(enc.ortho, mp, by='pahari_id') %>%
	filter(Chromosome!='X', Chromosome!='Y')


##hierarchical clusterting

expr = as.matrix(dplyr::select(enc.mp,
	cortex_encode_1, cortex_encode_2, Brain_pahari,
	liver_encode_1, liver_encode_2, Liver_pahari,
	kidney_encode_1, kidney_encode_2, Kidney_pahari,
	spleen_encode_1, spleen_encode_2, Spleen_pahari,
	placenta_encode_1, placenta_encode_2, Placenta_pahari_2, Placenta_pahari_3))
rownames(expr) = enc.mp$Probe

hc = hclust(dist(t(expr)))
par(mar=c(4,2,2,8))
plot(as.dendrogram(hc, hang=0.1), horiz=T)


##PCA on 5% most variable genes (~1000)

expr.var = apply(expr, 1, var)

species = case_when(grepl('encode', colnames(expr)) ~ 16,
					grepl('pahari', colnames(expr)) ~ 17)
										
tissue = rep(c('#648FFF','#785EF0','#DC267F','#FE6100','#FFB000'), c(rep(3,4),4))

pca = prcomp(expr[expr.var>quantile(expr.var,0.95),])

ggplot(pca$rotation, aes(x=PC1, y=PC2)) +
	geom_point(colour=tissue, shape=species, size=2) +
	theme_classic()

ggplot(pca$rotation, aes(x=PC2, y=PC3)) +
	geom_point(colour=tissue, shape=species, size=2) +
	theme_classic()



##expression of genes proximal to TEs active in the placenta

expr.av = tibble(gene = enc.mp$Probe,
				brain_musculus = rowMeans(select(enc.mp, cortex_encode_1:cortex_encode_2)),
				kidney_musculus = rowMeans(select(enc.mp, kidney_encode_1:kidney_encode_2)),
				liver_musculus = rowMeans(select(enc.mp, liver_encode_1:liver_encode_2)),
				placenta_musculus = rowMeans(select(enc.mp, placenta_encode_1:placenta_encode_2)),
				spleen_musculus = rowMeans(select(enc.mp, spleen_encode_1:spleen_encode_2)),
				brain_pahari = enc.mp$Brain_pahari,
				kidney_pahari = enc.mp$Kidney_pahari,
				liver_pahari = enc.mp$Liver_pahari,
				placenta_pahari = rowMeans(select(enc.mp, Placenta_pahari_2:Placenta_pahari_3)),
				spleen_pahari = enc.mp$Spleen_pahari) %>%
	pivot_longer(cols=!gene, names_to='sample', values_to='expr') %>%
	mutate(tissue = sub('_[[:print:]]+','',sample), species = sub('[[:print:]]+_','',sample)) %>%
	select(!sample) %>% pivot_wider(names_from = species, values_from = expr)

expr.diff = filter(expr.av, musculus > -1 | pahari > -1) %>%
	mutate(log2FC = musculus-pahari)


plot.te.genes = function(data, gene_file) {
	gene.data = read_tsv(gene_file)
	sub = filter(data, gene %in% gene.data$gene)
	p =	ggplot(sub, aes(x=tissue, y=log2FC)) +
			geom_boxplot(fill='lightblue', outlier.size=0.5) +
			theme_classic() +
			xlab('') +
			ylab('Log2 musculus/pahari') +
			scale_x_discrete(guide = guide_axis(angle = 90)) +
			geom_hline(yintercept=0, linetype='dashed')
	print(p)

	pval = group_by(sub, tissue) %>% summarise(p = t.test(log2FC)$p.value) %>%
		mutate(padj = p.adjust(p, method='BH'))
	return(pval)
}	


plot.te.genes(expr.diff,'IAP_proximal_genes.txt')
plot.te.genes(expr.diff,'RMER10A_mm_proximal_genes.txt')
plot.te.genes(expr.diff,'RMER10A_mp_proximal_genes.txt')





