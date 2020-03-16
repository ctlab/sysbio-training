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


## Using Combat to remove batch effect
batch <- pData(gse129260)$Replicate
modcombat <- model.matrix(~1, data=pData(gse129260))
combat_gse129260 = ComBat(dat=exprs(gse129260), batch=batch, mod=modcombat)


## PCA plot after batch removal
pcas <- prcomp(t(combat_gse129260), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(gse129260))
ggplot(plotData, aes(x=PC1, y=PC2, color=Cell, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)


## DE with the good model

cell_full_model <- model.matrix(~0 + Cell + Treatment + Replicate, data=pData(gse129260))
colnames(cell_full_model) <- c("il10neg", "il10pos", "LPS", "rep2")

fit <- lmFit(gse129260, cell_full_model)

fit2 <- contrasts.fit(fit, makeContrasts(il10pos - il10neg, levels=cell_full_model))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")


## Fancy volcano
ggplot(de, aes(x=logFC, y=-log10(adj.P.Val), color=adj.P.Val < 0.05)) +
  geom_point() + theme_bw() + scale_color_manual(values=c("black", "red")) +
  geom_text_repel(data=de %>% dplyr::filter(adj.P.Val < 0.01), aes(label=Gene.symbol, color=NULL))



## Performing fgsea ONLY
## kegg is the database with curated metabolic pathways

load("keggSymbolMouse.rdata")
stats <- de$t
names(stats) <- de$Gene.symbol


## Running fgsea
fgseaResults <- fgseaMultilevel(keggSymbolMouse, stats, minSize = 15, maxSize = 500)
head(fgseaResults, 3)

## Getting best pathways
topPathwaysUp <- fgseaResults[ES > 0, ][head(order(pval), n=5), pathway]
topPathwaysDown <- fgseaResults[ES < 0, ][head(order(pval), n=5), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))


## Plotting table for top pathways
dev.off()
plotGseaTable(keggSymbolMouse[topPathways], stats, fgseaResults, gseaParam = 0.5)

