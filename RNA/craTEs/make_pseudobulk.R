library(Matrix)
library(SingleCellExperiment)
library(glmGamPoi)


##Mouse gastrulation scRNA-seq data from:
##https://github.com/MarioniLab/EmbryoTimecourse2018

counts = readMM("raw_counts.mtx")
genes = read.table("genes.tsv", stringsAsFactors = F)
meta = read.table("meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")

rownames(counts) = genes[,1] #ensembl
colnames(counts) = meta$cell

sce = SingleCellExperiment(assays = list("counts" = counts))

colData(sce) = DataFrame(sample = meta$sample,
						stage = meta$stage,
						celltype = meta$celltype,
						row.names = meta$cell)


exe.epi = sce[,sce$celltype %in% c('Epiblast', 'ExE ectoderm')]

pb = pseudobulk(exe.epi, group_by = vars(celltype))

pb.mat = data.frame(counts(pb))

write.table(pb.mat, 'epi_exe_pseudobulk.txt', sep='\t', quote=F, col.names=NA)


##The Marsh et al pseudobulk was made at the end of 
##the script analysing those data in RNA/scRNA_processing
