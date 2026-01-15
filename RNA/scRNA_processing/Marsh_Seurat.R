library(Seurat)


##Read data from cellranger

plac.data = Read10X('filtered_feature_bc_matrix')


##Make and filter Seurat object

plac = CreateSeuratObject(plac.data, min.cells=3, min.features=200)

VlnPlot(plac, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plac = subset(plac, subset = nFeature_RNA > 200 & nFeature_RNA < 2000)


##Process

plac = NormalizeData(plac) 
plac = FindVariableFeatures(plac, selection.method = "vst", nfeatures = 2000)
plac = ScaleData(plac)
plac = RunPCA(plac)
ElbowPlot(plac, ndims=30)

plac = FindNeighbors(plac, dims = 1:20)
plac = FindClusters(plac, resolution = 0.5)


##Remove minor clusters

ncells = table(Idents(plac))
plac = subset(plac, idents = names(ncells)[ncells>50])


##Identify clusters

markers = FindAllMarkers(plac, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

VlnPlot(plac, features = c('Syna','Glis1','Stra6'), pt.size=0.5) #SynTI
VlnPlot(plac, features = c('Gcm1','Synb','Gcgr'), pt.size=0) #SynTII
VlnPlot(plac, features = c('Nos1ap','Lepr','Ctsq'), pt.size=0) #S-TGC
VlnPlot(plac, features = c('Prune2','Ncam1','Plac8'), pt.size=0) #Glycogen cells
VlnPlot(plac, features = c('Mitf','Flt1','Prl8a9'), pt.size=0) #SpT
VlnPlot(plac, features = c('Pecam1','Kdr','Tek'), pt.size=0.5) #Endothelial
VlnPlot(plac, features = c('Ebf1','Pdgfrb','Notch3'), pt.size=0) #Pericytes
VlnPlot(plac, features = c('Pdgfra','Kit'), pt.size=0) #Mesenchymal I
VlnPlot(plac, features = c('Col1a1','Col1a2'), pt.size=0) #Mesenchymal II
VlnPlot(plac, features = c('Pbx1','Pgr','Xist'), pt.size=0) #Decidua

new.cluster.ids = c("SynTI", "Glycogen cells", "SpT", "S-TGC", "Endothelial",
	"SynTII", "Pericytes", "Mesenchymal I", "Mesenchymal II", "Decidua")
names(new.cluster.ids) = levels(plac)
plac = RenameIdents(plac, new.cluster.ids)


##UMAP

plac = RunUMAP(plac, dims = 1:20)
DimPlot(plac, reduction = "umap", label=TRUE, repel=TRUE)


##Save Seurat object

SaveSeuratRds(plac, 'Marsh_IAP_Seurat.Rds')


##Make pseudobulk for craTEs analysis (in RNA/craTEs/)

pseudo = AggregateExpression(plac, assays = 'RNA')
write.table(as.matrix(pseudo$RNA), file='Marsh_E14_pseudobulk.txt', sep='\t', quote=F, col.names=NA)


