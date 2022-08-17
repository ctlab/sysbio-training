## loadin libraries and seurat object

library(Seurat)
library(Matrix)
library(MAST)
library(ggplot2)
library(dplyr)
library(fgsea)

seurat <- readRDS("blood_seurat.rds")



## calculate averaged expression per cluster
average <- AverageExpression(seurat)$SCT
averageLog <- log2(as.matrix(average) + 1)
colnames(averageLog) <- paste0("Cluster ", colnames(average))
write.table(averageLog, "average_log.tsv", sep="\t", col.names=NA, quote=F)


# Showing expression of some T cell populations markers

FeaturePlot(seurat, features=c("CD3D", "CD4", "CD8A", "PRF1"), ncol = 4)

VlnPlot(seurat, features=c("CD3D", "CD4", "CD8A", "PRF1"), ncol = 4, pt.size = 0.02)


## Differential expression cluster 3 vs cluster 11

de_03_vs_11 <- FindMarkers(
  seurat, assay="SCT", ident.1 = 3, ident.2 = 11,
  test="MAST", logfc.threshold = 0, min.pct = 0
)
write.table(de_03_vs_11, "de_03_vs_11.tsv", sep="\t", col.names=NA, quote=F)
topGenes <- head(rownames(de_03_vs_11))



## showing how table looks like
head(de_03_vs_11)


## violin plot for some top genes of DE
VlnPlot(seurat, topGenes, pt.size = 0.02, idents=c(3, 11), ncol=6)


## saving top genes of DE into text files for further pathway enrichment analysis
de_03_vs_11$gene <- rownames(de_03_vs_11)

top50 <- de_03_vs_11 %>% top_n(50, avg_log2FC) %>% pull(gene)
top200 <- de_03_vs_11 %>% top_n(200, avg_log2FC) %>% pull(gene)
bottom50 <- de_03_vs_11 %>% top_n(50, -avg_log2FC) %>% pull(gene)
bottom200 <- de_03_vs_11 %>% top_n(200, -avg_log2FC) %>% pull(gene)

writeLines(top50, "top_50.txt")
writeLines(top200, "top_200.txt")
writeLines(bottom50, "bottom_50.txt")
writeLines(bottom200, "bottom_200.txt")


## fgsea pathway enrichment against kegg database
load("keggSymbolHuman.rdata")

ranks <- de_03_vs_11$avg_log2FC
names(ranks) <- rownames(de_03_vs_11)
fgseaRes <- fgsea(pathways = keggSymbolHuman, 
                  stats = ranks,
                  minSize=15,
                  maxSize=500,
                  nperm=100000)


## showing how fgsea results look like
head(fgseaRes)



## choosing top pathways
topPathwaysUp <- fgseaRes[ES > 0 & padj < 0.01, ][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0 & padj < 0.01, ][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))


## plotting top pathways
plotGseaTable(keggSymbolHuman[topPathways], ranks, fgseaRes, 
              gseaParam = 0.2, colwidths = c(5, 1, 0.8, 0.8, 0.8))

## showing enrichment for one of the pathways
plotEnrichment(keggSymbolHuman[["T cell receptor signaling pathway - Homo sapiens (human)"]],
ranks) + labs(title="T cell receptor signaling pathway - Homo sapiens (human)")

