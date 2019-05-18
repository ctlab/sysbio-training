library(GEOquery)
library(limma)
library(org.Mm.eg.db)

# for collapseBy
source("functions.R")

dir.create("cache")

es <- getGEO("GSE53986", AnnotGPL = TRUE, destdir="cache")[[1]]
str(experimentData(es))
str(pData(es))
head(fData(es))
es$`treatment:ch1`

es$condition <- gsub("\\+", "_", es$`treatment:ch1`)
es$condition

es <- collapseBy(es, fData(es)$`Gene ID`, FUN=median)
es <- es[!grepl("///", rownames(es)), ]
es <- es[rownames(es) != "", ]

# there is a lot of garbage there
fData(es) <- data.frame(row.names = rownames(es))
fData(es)$symbol <- mapIds(org.Mm.eg.db, keys = rownames(es), column = "SYMBOL", keytype = "ENTREZID")

es.qnorm <- es

summary(exprs(es.qnorm))
exprs(es.qnorm) <- normalizeBetweenArrays(log2(exprs(es.qnorm)+1), method="quantile")
summary(exprs(es.qnorm))

es.qnorm.top12K <- es.qnorm
es.qnorm.top12K <- es.qnorm.top12K[head(order(apply(exprs(es.qnorm.top12K), 1, mean), 
                                              decreasing = TRUE), 12000), ]

es.design <- model.matrix(~0+condition, data=pData(es.qnorm.top12K))

fit <- lmFit(es.qnorm.top12K, es.design)

fit2 <- contrasts.fit(fit, makeContrasts(conditionLPS-conditionUntreated, 
                                         levels=es.design))

# fit2 <- contrasts.fit(fit, makeContrasts2(c("condition", "LPS", "Untreated"),
#                                          levels=es.design))
fit2 <- eBayes(fit2)
de <- topTable(fit2, adjust.method="BH", number=Inf)
head(de)
library(data.table)
de <- as.data.table(de, keep.rownames=TRUE)
de[symbol == "Nos2"]
