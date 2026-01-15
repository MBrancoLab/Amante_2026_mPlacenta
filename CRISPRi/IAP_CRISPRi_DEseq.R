library('DESeq2')
library('tidyverse')
library('topGO')
library('org.Mm.eg.db')


##read raw counts

counts = read.delim('gene_raw_counts.txt')

reads = as.matrix(counts[,13:ncol(counts)])
rownames(reads) = counts$Probe
colData = DataFrame(sample=factor(rep(c('KRAB','IAP_gRNA'),each=3)))


##get DE genes

dds = DESeqDataSetFromMatrix(reads,colData,design=~sample)
dds = DESeq(dds)
res = results(dds,contrast=c('sample','IAP_gRNA','KRAB'))
up = subset(res, padj<0.05 & log2FoldChange>1)
down = subset(res, padj<0.05 & log2FoldChange<(-1))


##Volcano plot

quartz(w=3,h=3)
mutate(as.data.frame(res), is.sig = padj<0.05 & abs(log2FoldChange)>1) %>%
ggplot(aes(x=log2FoldChange, y=-log10(padj), colour=is.sig)) +
	geom_point(size=0.5) +
	theme_classic() +
	scale_colour_manual(values=c('grey','blue'))


##gene ontology analysis

go.ann = function(de.genes) {
	genes = factor(as.integer(rownames(reads) %in% rownames(de.genes)))
	names(genes) = rownames(reads)
	go = new('topGOdata', ontology = 'BP',
		allGenes = genes, geneSelectionFun = function(x) x==0,
		annotationFun = annFUN.org, mapping = 'org.Mm.eg.db', ID = 'SYMBOL')
	return(go)
}

go.plot = function(go, nodes=8) {
	rt = runTest(go, algorithm='classic', statistic='fisher')
	gt = GenTable(go, classicFisher=rt, topNodes=nodes)

	ggplot(gt, aes(x=-log10(as.numeric(classicFisher)), y=Term)) +
		geom_bar(stat='identity', fill='grey', colour='black', width=0.7) +
		theme_classic() +
		ylab('') +
		xlab('-log10 p')
}

go.up = go.ann(up)
go.plot(go.up)

go.down = go.ann(down)
go.plot(go.down)

