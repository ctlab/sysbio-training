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

# creating an expression set object

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

# PCA

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


p <- pcaPlot(es.qnorm.top12K, 1, 2) + 
  aes(color=condition) + 
  geom_text_repel(aes(label=sample)) 
print(p)
ggsave(p, width=6, height=4, file="./es.pca12.png")

# differential expression with DESeq2

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
# log2FC = IL4 - untr
de <- results(dds, contrast = c("condition", "IL4", "untr"), cooksCutoff = F)
head(de)
de <- data.table(ID=rownames(de), as.data.table(de))
head(de)

head(fData(es))
de <- cbind(de, fData(es))
de <- de[ID %in% rownames(es.qnorm.top12K), ]
de <- de[order(stat), ]
de

de[symbol == "Arg1"]

fwrite(de, file="./de/untr.vs.IL4.de.tsv", sep="\t")


# pathway analysis with fgsea

stats <- de[, setNames(stat, entrez)]
library(msigdbr)
# GO BP pathways from MSigDB
m_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")
m_df
pathways <- split(m_df$entrez_gene, m_df$gs_name)

library(fgsea)
fr <- fgsea(pathways, stats, nperm = 100000, nproc=4, minSize=15, maxSize=500)
fr[order(pval)]
frML <- fgseaMultilevel(pathways, stats,
                        sampleSize = 100,
                        nproc=4, minSize=15, maxSize=500)
frML[order(pval)]

frML[padj < 0.01]

collapsedPathways <- collapsePathways(fr[order(pval)][padj < 0.01], pathways, stats)
str(collapsedPathways)

mainPathways <- frML[pathway %in% collapsedPathways$mainPathways][
  order(sign(ES)*log(pval)), pathway]

frMain <- frML[match(mainPathways, pathway)]
frMain[, leadingEdge := lapply(leadingEdge, mapIds, 
                               x=org.Mm.eg.db, keytype="ENTREZID", column="SYMBOL")]
dir.create("gsea")
fwrite(frMain, file="gsea/untr.vs.IL4.filtered.tsv", sep="\t", sep2=c("", " ", ""))

pdf("gsea/untr.vs.IL4.pdf", width=12, height=2 + length(mainPathways) * 0.25)
plotGseaTable(pathways = pathways[mainPathways], stats = stats, fgseaRes=frMain, gseaParam = 0.5)
dev.off()

fr[, leadingEdge := NULL]
fwrite(frML[order(sign(ES)*log(pval))], file="./gsea/untr.vs.IL4.full.tsv", sep="\t") 

plotEnrichment(pathways[["GO_DEFENSE_RESPONSE_TO_VIRUS"]], stats) + 
  ggtitle("Defense reponse to virus")

plotEnrichment(pathways[["GO_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS"]], stats) + 
  ggtitle("Nucleotide triphosphate metabolism")

plotEnrichment(pathways[["GO_CELLULAR_RESPONSE_TO_CYTOKINE_STIMULUS"]], stats) + 
  ggtitle("Cellular response to cytokines")
