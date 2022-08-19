## loading libraries
library(Seurat)
library(Matrix)
library(MAST)
library(ggplot2)
library(dplyr)


## checking dataset dimensions
data <- Read10X("/home/student/shared/filtered_feature_bc_matrix/")
dim(data)


## looking at how UMIs are distributed

plotData <- data.frame(
  umis <- colSums(data)
)
ggplot(data=plotData, aes(x=umis)) +
  geom_histogram() + theme_bw()



## creating seurat object
seurat <- CreateSeuratObject(data, min.cells = 10, min.features = 10)
dim(seurat)



## UMIs vs Genes detected
FeatureScatter(seurat, "nCount_RNA", "nFeature_RNA") + scale_x_log10() + scale_y_log10()


## mitochondrial content
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
FeatureScatter(seurat, "nCount_RNA", "percent.mt") + scale_x_log10()
FeatureScatter(seurat, "nFeature_RNA", "percent.mt") + scale_x_log10()


## filtering bad cells
seurat <- subset(seurat, subset = nFeature_RNA > 1000 & percent.mt < 15)
dim(seurat)


## old way to normalize data
## seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)
## seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
## seurat <- ScaleData(seurat)



## normalize in one command
seurat <- SCTransform(seurat, vars.to.regress = "percent.mt", verbose = FALSE)



## Plot showing variable features
VariableFeaturePlot(seurat) + scale_y_log10()


## Plot showing variable features
top10_variable_genes <- head(VariableFeatures(seurat), 10)
VariableFeaturePlot(seurat) %>% 
  LabelPoints(points = top10_variable_genes, repel = TRUE) +
  scale_y_log10()


## PCA and elbow plot
seurat <- RunPCA(seurat, verbose = FALSE)
ElbowPlot(seurat, ndims = 50)


## TSNE
seurat <- RunTSNE(seurat, dims=1:20)
DimPlot(seurat, reduction = "tsne") + NoLegend()


## UMAP
seurat <- RunUMAP(seurat, dims=1:20)
DimPlot(seurat, reduction = "umap") + NoLegend()


## comparing TSNE vs UMAP
DimPlot(seurat, reduction = "tsne") + NoLegend()
DimPlot(seurat, reduction = "umap") + NoLegend()


## clustering
seurat <- FindNeighbors(seurat, dims = 1:20, verbose = FALSE)
seurat <- FindClusters(seurat, resolution=0.6, verbose = FALSE)
DimPlot(seurat, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()


## showing gene expression on top of dimensionality reduction plot
FeaturePlot(seurat, c("CD14", "CD79A", "CD3D"), cols=c("grey", "red"), reduction="umap", ncol=3)


## finding markers for all clusters
## max cells per ident is only seed to speed up the whole thing
allMarkers <- FindAllMarkers(seurat, max.cells.per.ident = 100, test.use = "MAST", only.pos = T)
goodMarkers <- allMarkers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC) %>% pull(gene)
goodMarkers


## showing gene expression on top of dimensionality reduction plot
FeaturePlot(seurat, goodMarkers[1:3], cols=c("grey", "red"), reduction="umap", ncol=3)


## showing gene expression as a violin plot
VlnPlot(seurat, goodMarkers[1:3], pt.size = 0.1)
