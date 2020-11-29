# Step 0: importing gene count matrix

source("./functions.R")
library(data.table)

# load the featureCount files
fc.files <- list.files("./featureCounts/", pattern="fc.txt$", recursive = TRUE, full.names = TRUE)
fc_res <- lapply(fc.files, fread)

# compaile to a matrix
fc_mat <- do.call(cbind, lapply(fc_res, function(x) x[[ncol(x)]]))
rownames(fc_mat) <- fc_res[[1]][, Geneid]
tags <- sapply(fc.files, function (x) gsub(".fc.txt", "", basename(x), fixed=TRUE))
colnames(fc_mat) <- tags

head(fc_mat)
tail(fc_mat)

# Step 1: creating an ExpressionSet object

library(Biobase)

es <- ExpressionSet(fc_mat)
pData(es)$condition <- c(rep("Ctrl", 3), rep("LPS", 3))

head(pData(es))
head(fData(es))

head(rownames(es))
fData(es)$ensembl <- gsub("\\.\\d*$", "", rownames(es))

library(org.Mm.eg.db)
columns(org.Mm.eg.db)
fData(es)$entrez <- mapIds(org.Mm.eg.db, keys=fData(es)$ensembl, 
                           keytype="ENSEMBL", column="ENTREZID" )
fData(es)$symbol <- mapIds(org.Mm.eg.db, keys=fData(es)$ensembl, 
                           keytype="ENSEMBL", column="SYMBOL" )

head(fData(es))

exprs(es)[which(fData(es)$symbol == "Actb"), ]
exprs(es)[which(fData(es)$symbol == "Acod1"), ]

# Step 2: getting a quantile normalized file

library(limma)
library(ggplot2)
library(ggrepel)

es.qnorm <- es

exprs(es.qnorm) <- normalizeBetweenArrays(log2(exprs(es.qnorm) + 1), method="quantile")

es.qnorm.top12K <- es.qnorm
fData(es.qnorm.top12K)$mean <- apply(exprs(es.qnorm.top12K), 1, mean)
es.qnorm.top12K <- es.qnorm.top12K[order(fData(es.qnorm.top12K)$mean, decreasing = TRUE), ]
head(exprs(es.qnorm.top12K))
es.qnorm.top12K <- es.qnorm.top12K[!duplicated(fData(es.qnorm.top12K)$entrez), ]
es.qnorm.top12K <- es.qnorm.top12K[!is.na(fData(es.qnorm.top12K)$entrez), ]
es.qnorm.top12K <- es.qnorm.top12K[1:12000,]
write.gct(es.qnorm.top12K, file="./es.qnorm.top12k.gct")

# Step 3: PCA plot

p <- pcaPlot(es.qnorm.top12K, 1, 2) + 
  aes(color=condition) + 
  geom_text_repel(aes(label=sample)) 
print(p)
ggsave(p, width=6, height=4, file="./es.pca12.png")

# Step 4: differential expression with DESeq2

library(DESeq2)
dds <- DESeqDataSetFromMatrix(exprs(es), pData(es), design=~condition)
dds
dds <- DESeq(dds)
dds

plotDispEsts(dds)
vst <- varianceStabilizingTransformation(dds)
plotPCA(vst)

dir.create("./de/", showWarnings = F)
unique(dds$condition)
# log2FC = LPS - Ctrl
de <- results(dds, contrast = c("condition", "LPS", "Ctrl"), cooksCutoff = F)
head(de)
de <- data.table(ID=rownames(de), as.data.table(de))
head(de)

head(fData(es))
de <- cbind(de, fData(es))
de <- de[ID %in% rownames(es.qnorm.top12K), ]
de <- de[order(stat), ]
de

de[symbol == "Acod1"]

fwrite(de, file="./de/Ctrl.vs.LPS.de.tsv", sep="\t")

# Step 5: pathway analysis with fgsea

stats <- de[, setNames(stat, entrez)]
library(msigdbr)
# Hallmark pathways from MSigDB
m_df <- msigdbr(species = "Mus musculus", category = "H")
m_df
pathways <- split(m_df$entrez_gene, m_df$gs_name)

library(fgsea)

fr <- fgseaMultilevel(pathways, stats,
                      sampleSize = 100,
                      nproc=4, 
                      minSize=15, maxSize=500, 
                      absEps=1e-10)
fr[padj < 0.01][order(NES)]

dir.create("gsea")
fwrite(fr, file="gsea/Ctrl.vs.LPS.filtered.tsv", sep="\t", sep2=c("", " ", ""))

mainPathways <- fr[padj < 0.01][order(NES)][, pathway]
pdf("gsea/Ctrl.vs.LPS.pdf", width=12, height=2 + length(mainPathways) * 0.25)
plotGseaTable(pathways = pathways[mainPathways], stats = stats, fgseaRes=fr, gseaParam = 0.5)
dev.off()


plotEnrichment(pathways$HALLMARK_INTERFERON_GAMMA_RESPONSE, stats = stats) +
  ggtitle("INFg response")
