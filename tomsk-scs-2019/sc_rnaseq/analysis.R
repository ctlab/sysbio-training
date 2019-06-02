library(Seurat)
library(ggplot2)
library(fgsea)
library(dplyr)

set.seed(1)
setwd("~/scRNAseq")

# PART 1: reading data and creating dataset object
counts <- Read10X("filtered_feature_bc_matrix")
pbmc <- CreateSeuratObject(counts = counts, project = "PBMC")

# umis <- data.frame(nUmi = pbmc@meta.data$nCount_RNA)
# ggplot(data=umis, aes(x=nUmi)) +
#   geom_histogram() + theme_bw(base_size = 8) + 
#   labs(y="Cell count", x="UMIs in the cell")
# ggsave("umis_total.pdf", width=3, height=2)


# PART 2: calculating mitochondrial content and quality checks
mito.genes <- c(grep("^MT-", rownames(x = pbmc), value = T),
                grep("^mt-", rownames(x = pbmc), value = T))

percent.mito <- Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts')[mito.genes, ]) / 
  Matrix::colSums(x = GetAssayData(object = pbmc, slot = 'counts'))
pbmc[['percent.mito']] <- percent.mito
pbmc[['nCount_RNA_log10']] <- log10(pbmc[['nCount_RNA']])
pbmc[['nFeature_RNA_log10']] <- log10(pbmc[['nFeature_RNA']])

FeatureScatter(object = pbmc, feature1 = "nCount_RNA_log10", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nCount_RNA_log10", feature2 = "nFeature_RNA_log10")

# PART 3: Filtering out unwanted cells and scaling dataset
pbmc <- subset(x = pbmc, subset = percent.mito < 0.25 & nFeature_RNA_log10 > 2.5)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

VariableFeaturePlot(pbmc)

# PART 4: Scaling data, takes some time
pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("nCount_RNA"))

# PART 4.1 genes detected within a cell
cell1 <- as.numeric(GetAssayData(object = pbmc, slot = 'data')[, 1])
umiData <- data.frame(gene=rownames(GetAssayData(object = pbmc, slot = 'data')), expression=cell1)
ggplot(data=umiData) +
  geom_histogram(aes(x=expression)) +
  scale_y_log10(breaks=c(1, 10, 100, 1000, 10000)) + theme_bw(base_size=8) +
  labs(y="Frequency", x="log2(UMIs + 1)")
ggsave("umis.pdf", width=3, height=2)


# PART 5: PCA and elbow plot

pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc), npcs = 50,
               ndims.print = 1:5, nfeatures.print = 10)
ElbowPlot(object = pbmc, ndims = 50)

# PART 5: dimensionality reductions. tSNE and uMAP
pbmc <- RunTSNE(pbmc, dims = 1:10)
pbmc <- RunUMAP(pbmc, dims = 1:10)


DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "umap")

FeaturePlot(pbmc, features = c("CD14", "CD79A", "CD3E", "CD8A"), reduction = "tsne")
FeaturePlot(pbmc, features = c("CD14", "CD79A", "CD3E", "CD8A"), reduction = "umap")


## PART 6: Clustering
pbmc <- FindNeighbors(object = pbmc, dims = 1:10)
pbmc <- FindClusters(object = pbmc, resolution = 0.6)

DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "umap")


## PART 7: Averaging expression per cluster

average <- AverageExpression(pbmc)$RNA
averageLog <- log2(as.matrix(average) + 1)
colnames(averageLog) <- paste0("Cluster ", colnames(average))
write.table(averageLog, "average_log.tsv", sep="\t", col.names=NA, quote=F)


# PART 8: comparing T cells.

FeaturePlot(pbmc, features = c("CD3E", "CD4", "CD8A", "PRF1"), reduction = "tsne")
VlnPlot(pbmc, c("CD3E", "CD4", "PRF1", "CD8A"), ncol = 1, pt.size = 0.02)

de_0_vs_3 <- FindMarkers(pbmc, ident.1 = 0, ident.2 = 3, test="MAST", latent.vars = "nCount_RNA", logfc.threshold = 0, min.pct = 0)
head(de_0_vs_3)

topGenes <- rownames(de_0_vs_3[1:6, ])
VlnPlot(pbmc, topGenes, pt.size = 0.02)
VlnPlot(pbmc, topGenes, pt.size = 0.02, idents = c(0, 3))
write.table(de_0_vs_3, "de_0_vs_3.tsv", sep="\t", col.names=NA, quote=F)

# PART 9: saving top genes

de_0_vs_3$gene <- rownames(de_0_vs_3)

top50 <- de_0_vs_3 %>% top_n(50, avg_logFC) %>% pull(gene)
top200 <- de_0_vs_3 %>% top_n(200, avg_logFC) %>% pull(gene)
bottom50 <- de_0_vs_3 %>% top_n(50, -avg_logFC) %>% pull(gene)
bottom200 <- de_0_vs_3 %>% top_n(200, -avg_logFC) %>% pull(gene)

writeLines(top50, "top_50.txt")
writeLines(top200, "top_200.txt")
writeLines(bottom50, "bottom_50.txt")
writeLines(bottom200, "bottom_200.txt")

# PART 10: pathway enrichment

load("~/scRNAseq/scexplorer/keggSymbolHuman.rdata")

ranks <- de_0_vs_3$avg_logFC
names(ranks) <- rownames(de_0_vs_3)
fgseaRes <- fgsea(pathways = keggSymbolHuman,
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)
head(fgseaRes)

graphics.off()
topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.01, ][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.01, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(keggSymbolHuman[topPathways], ranks, fgseaRes,
              gseaParam = 0.2, colwidths = c(5, 1, 0.8, 0.8, 0.8))


plotEnrichment(keggSymbolHuman[["T cell receptor signaling pathway - Homo sapiens (human)"]],
               ranks) + labs(title="T cell receptor signaling pathway - Homo sapiens (human)")

## PART 11: Finding cluster markers
dir.create("scexplorer/data/10x_1k", showWarnings = F)
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.2, thresh.use = 0.10, test.use = "MAST", latent.vars = "nCount_RNA")
write.table(pbmc.markers, "scexplorer/data/10x_1k/markers.tsv", sep="\t", quote=F, row.names=F)
 
# PART 12: saving data for explorer

## SAVING DATA FOR EXPLORER

expData <- GetAssayData(object = pbmc, slot = 'data')
save(expData, file="scexplorer/data/10x_1k/expData.Rda")

dataForPlot <- as.data.frame(pbmc@reductions$tsne@cell.embeddings)
dataForPlot$Cluster <-  Idents(object = pbmc)
dataForPlot$nUmi <- pbmc@meta.data$nCount_RNA
dataForPlot$nGene <- pbmc@meta.data$nFeature_RNA
dataForPlot$nUmiLog2 <- log2(pbmc@meta.data$nCount_RNA)
dataForPlot$nGeneLog2 <- log2(pbmc@meta.data$nFeature_RNA)

write.table(dataForPlot, "scexplorer/data/10x_1k/data_for_plot.tsv", sep="\t", quote=F, col.names=NA)
