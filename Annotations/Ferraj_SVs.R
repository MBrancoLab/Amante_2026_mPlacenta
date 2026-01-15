library(tidyverse)


##Ferraj et al supp data filtered for IAPs with H3K27ac in placenta

del = read_tsv('active_IAP_dels_mm39.txt')



##parse strain names into matrix of deletion calls

samples = strsplit(del$MERGE_SAMPLES, split=',')
strains = unique(unlist(samples))

isdel = lapply(samples, function(x) strains %in% x)
delmat = matrix(as.numeric(unlist(isdel)), nrow=length(isdel), byrow=TRUE)
colnames(delmat) = strains
rownames(delmat) = del$ID


##plot as heatmap

h = heatmap(delmat, scale='none', col=c('navy','white'))