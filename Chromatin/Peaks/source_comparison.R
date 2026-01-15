library(gplots)


##Count matrix (from Seqmonk, using merged ENCODE peaks)

counts = read.delim('H3K27ac_log2Counts.txt.gz')


##Make correlation matrix

nc = ncol(counts)
cor.mat = matrix(nrow=nc-12, ncol=nc-12)
colnames(cor.mat) = rownames(cor.mat) = colnames(counts)[13:nc]
for (i in 13:nc) {
	for (j in 13:nc) {
		cor.mat[i-12,j-12] = cor.test(counts[,i],counts[,j])$estimate
	}
}


##Plot

heatmap.2(cor.mat, scale='none',
	col = colorRampPalette(c('navy','black','yellow')),
	breaks = seq(0.5,0.85,length.out=100),
	margins = c(8,8), cexRow=0.7, cexCol=0.7,
	trace = 'none', density.info='none')
