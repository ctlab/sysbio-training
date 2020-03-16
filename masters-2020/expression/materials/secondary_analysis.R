## Installing packages if needed

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("sva", quietly = TRUE)) BiocManager::install("sva")
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")


## loading packages
library(Biobase)
library(limma)
library(sva)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(fgsea)

# setting up color pallete
# loading data

blueWhiteRed <- colorRampPalette(c("#3859A8", "#EEEEEE", "#EE2930"))(10)

load("gse129260.Rdata")




## Looking only at 3 genes for now

someGenes <- exprs(gse129260)[c("Actb", "Ddx3y", "Il10"), ]
plotData <- t(someGenes)
plotData <- as.data.frame(plotData)
plotData <- cbind(plotData, pData(gse129260))

head(plotData, 4)



## All the plots for these 3 genes are not neccessary analysis steps
## Plots: Actb
ggplot(plotData, aes(y=Actb)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Cell, y=Actb)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Replicate, y=Actb)) +
  geom_boxplot() + theme_bw()



## Plots: Il10
ggplot(plotData, aes(y=Il10)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Cell, y=Il10)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Replicate, y=Il10)) +
  geom_boxplot() + theme_bw()


## Plots: Ddx3y
ggplot(plotData, aes(y=Ddx3y)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Cell, y=Ddx3y)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Replicate, y=Ddx3y)) +
  geom_boxplot() + theme_bw()


## Heatmap for those genes
pheatmap(someGenes, scale="row", color=blueWhiteRed, annotation_col = pData(gse129260), cluster_cols = F)



## Clustered heatmap
pheatmap(exprs(gse129260), scale="row", color=blueWhiteRed, border_color = NA, kmeans_k = 8,
         annotation_col = pData(gse129260), cluster_cols = F)



## Running PCA
pcas <- prcomp(t(exprs(gse129260)), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(gse129260))


## Plotting PCA
ggplot(plotData, aes(x=PC1, y=PC2, color=Cell, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)
ggplot(plotData, aes(x=PC1, y=PC2, color=Replicate, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)


## Identifying genes "parallel" to PC1
rotation <- pcas$rotation
PC1GenesDown <- head(rownames(rotation[order(rotation[, 1]), ]), 10)
PC1GenesUp <- tail(rownames(rotation[order(rotation[, 1]), ]), 10)
print(PC1GenesDown)
print(PC1GenesUp)


## Showing heatmap for those genes
pheatmap(exprs(gse129260)[c(PC1GenesDown, PC1GenesUp), ], 
         scale="row", color=blueWhiteRed, border_color = NA,
         annotation_col = pData(gse129260), cluster_cols = F)



## Identifying genes "parallel" to PC2
rotation <- pcas$rotation
PC2GenesDown <- head(rownames(rotation[order(rotation[, 2]), ]), 10)
PC2GenesUp <- tail(rownames(rotation[order(rotation[, 2]), ]), 10)
print(PC2GenesDown)
print(PC2GenesUp)


## Showing heatmap for those genes
pheatmap(exprs(gse129260)[c(PC2GenesDown, PC2GenesUp), ], 
         scale="row", color=blueWhiteRed, border_color = NA,
         annotation_col = pData(gse129260), cluster_cols = F)



## Using Combat to remove batch effect
batch <- pData(gse129260)$Replicate
modcombat <- model.matrix(~1, data=pData(gse129260))
combat_gse129260 = ComBat(dat=exprs(gse129260), batch=batch, mod=modcombat)


## PCA plot after batch removal
pcas <- prcomp(t(combat_gse129260), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(gse129260))
ggplot(plotData, aes(x=PC1, y=PC2, color=Cell, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)


## Same 3 genes as before, but without batch

someGenes <- combat_gse129260[c("Actb", "Ddx3y", "Il10"), ]
plotData <- t(someGenes)
plotData <- as.data.frame(plotData)
plotData <- cbind(plotData, pData(gse129260))

head(plotData)




## Actb
ggplot(plotData, aes(y=Actb)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Cell, y=Actb)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Replicate, y=Actb)) +
  geom_boxplot() + theme_bw()



## Il10
ggplot(plotData, aes(y=Il10)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Cell, y=Il10)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Replicate, y=Il10)) +
  geom_boxplot() + theme_bw()


## Ddx3y
ggplot(plotData, aes(y=Ddx3y)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Cell, y=Ddx3y)) +
  geom_boxplot() + theme_bw()
ggplot(plotData, aes(x=Replicate, y=Ddx3y)) +
  geom_boxplot() + theme_bw()

## Below are models

## ------------------------------------------------------------------------
model_simple <- model.matrix(~0 + Cell, data = pData(gse129260))
colnames(model_simple) <- c("Negative", "Positive")
model_simple


## ------------------------------------------------------------------------
exprs(gse129260)["Il10", ]
linear_fit <- lm.fit(model_simple, exprs(gse129260)["Il10", ])
linear_fit$coef


## ----fig.show = 'hold', fig.width=6, fig.height=2, dev='svg'-------------
pheatmap(exprs(gse129260)[c("Il10", "S100a10"), ], scale="row", cluster_cols = F, cluster_rows = F, color=blueWhiteRed, annotation_col = pData(gse129260))


## ------------------------------------------------------------------------
treatment_model <- model.matrix(~0 + Treatment, data = pData(gse129260))
colnames(treatment_model) <- c("aCD-40", "LPS")
treatment_model


## ------------------------------------------------------------------------
treatment_rep_model <- model.matrix(~0 + Treatment + Replicate, data = pData(gse129260))
colnames(treatment_rep_model) <- c("aCD-40", "LPS", "Rep2")
treatment_rep_model


## ------------------------------------------------------------------------
exprs(gse129260)["S100a10", ]
linear_fit <- lm.fit(treatment_rep_model, exprs(gse129260)["S100a10", ])
linear_fit$coef


## ------------------------------------------------------------------------
full_model <- model.matrix(~0 + Treatment + Cell + Replicate, data = pData(gse129260))
colnames(full_model) <- c("aCD-40", "LPS", "Il10pos", "Rep2")
full_model


## ------------------------------------------------------------------------
linear_fit <- lm.fit(full_model, exprs(gse129260)["Il10", ])
linear_fit$coef

linear_fit <- lm.fit(full_model, exprs(gse129260)["S100a10", ])
linear_fit$coef


## ------------------------------------------------------------------------
full_model <- model.matrix(~1 + Treatment + Cell + Replicate, data = pData(gse129260))
colnames(full_model) <- c("Intercept", "LPS", "Il10pos", "Rep2")
full_model


## ------------------------------------------------------------------------
linear_fit <- lm.fit(full_model, exprs(gse129260)["Il10", ])
linear_fit$coef

linear_fit <- lm.fit(full_model, exprs(gse129260)["S100a10", ])
linear_fit$coef


## DE with the good model

cell_full_model <- model.matrix(~0 + Cell + Treatment + Replicate, data=pData(gse129260))
colnames(cell_full_model) <- c("il10neg", "il10pos", "LPS", "rep2")

fit <- lmFit(gse129260, cell_full_model)

fit2 <- contrasts.fit(fit, makeContrasts(il10pos - il10neg, levels=cell_full_model))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## ------------------------------------------------------------------------
head(de)


## Volcano plot for DE
ggplot(de, aes(x=logFC, y=-log10(adj.P.Val), color=adj.P.Val < 0.05)) +
  geom_point() + theme_bw()




## Fancy volcano
ggplot(de, aes(x=logFC, y=-log10(adj.P.Val), color=adj.P.Val < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red")) +
  geom_text_repel(data=de %>% dplyr::filter(adj.P.Val < 0.01), aes(label=Gene.symbol, color=NULL))




## DE with the bad model
cell_bad_model <- model.matrix(~0 + Cell, data=pData(gse129260))
colnames(cell_bad_model) <- c("il10neg", "il10pos")

fit_bad <- lmFit(gse129260, cell_bad_model)

fit_bad2 <- contrasts.fit(fit_bad, makeContrasts(il10pos - il10neg, levels=cell_bad_model))
fit_bad2 <- eBayes(fit_bad2, trend = T)

de_bad <- topTable(fit_bad2, adjust.method="BH", number=Inf, sort.by = "P")


## Comparing DE results (in # of genes)
de %>% filter(adj.P.Val < 0.05) %>% count()
de_bad %>% filter(adj.P.Val < 0.05) %>% count()



## Performing exact Fisher Test
## kegg is the database with curated metabolic pathways

load("keggSymbolMouse.rdata")
upRegulatedGenes <- de %>% filter(adj.P.Val < 0.05 & logFC > 0) %>% pull("Gene.symbol")
length(upRegulatedGenes)

randomGeneSet <- keggSymbolMouse[["Cardiac muscle contraction - Mus musculus (mouse)"]]
randomGeneSet <- randomGeneSet[randomGeneSet %in% rownames(de)]
length(randomGeneSet)

length(intersect(randomGeneSet, upRegulatedGenes))



## Performing the EXACT fisher test (for random unrelated gene set)

N <- nrow(de)
K <- length(randomGeneSet)
n <- length(upRegulatedGenes)
k <- length(intersect(upRegulatedGenes, randomGeneSet))
phyper(k - 1, K, N - K, n, lower.tail = F)


## Performing for non-random set

nonRandomGeneSet <- keggSymbolMouse[["Cytokine-cytokine receptor interaction - Mus musculus (mouse)"]]
nonRandomGeneSet <- nonRandomGeneSet[nonRandomGeneSet %in% rownames(de)]


N <- nrow(de)
K <- length(nonRandomGeneSet)
n <- length(upRegulatedGenes)
k <- length(intersect(upRegulatedGenes, nonRandomGeneSet))
print(c(N, K, n, k))
phyper(k - 1, K, N - K, n, lower.tail = F)


## Paste these gene in MsigDB
cat(upRegulatedGenes)


## Enrichment plot: random gene set
stats <- de$t
names(stats) <- de$Gene.symbol
plotEnrichment(randomGeneSet, stats)



## Enrichmet plot: non-random gene set
plotEnrichment(nonRandomGeneSet, stats)



## Running fgsea
fgseaResults <- fgseaMultilevel(keggSymbolMouse, stats, minSize = 15, maxSize = 500)
head(fgseaResults, 3)



## Getting best pathways
topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=5), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))


## Plotting table for top pathways
dev.off() ## might be usefull to remove previous plots
plotGseaTable(keggSymbolMouse[topPathways], stats, fgseaResults, gseaParam = 0.5)

