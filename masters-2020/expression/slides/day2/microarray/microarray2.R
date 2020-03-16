## ----message=F-----------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("sva", quietly = TRUE)) BiocManager::install("sva")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")




library(GEOquery)
library(Biobase)
library(ggplot2)
library(reshape2)
library(limma)
library(MASS)


## ----cache=T, message=F--------------------------------------------------

GSE129260 <- getGEO("GSE129260", AnnotGPL = TRUE)[[1]]

pData(GSE129260)


## ----message=F-----------------------------------------------------------
pData(GSE129260)$rep <- gsub(".*(rep\\d)$", "\\1", pData(GSE129260)$title)
pData(GSE129260) <- pData(GSE129260)[, c("characteristics_ch1.1", "characteristics_ch1.2", "rep")]

colnames(pData(GSE129260)) <- c("Cell", "Treatment", "Replicate")
head(pData(GSE129260))



## ----cache=T, message=F--------------------------------------------------
colnames(fData(GSE129260))



## ----message=F-----------------------------------------------------------
fData(GSE129260) <- fData(GSE129260)[, c("ID", "Gene symbol", "Gene ID")]
head(fData(GSE129260))



## ----fig.height=3, fig.fullwidth=T, dev='svg'----------------------------

ggplot(data=data.frame(expression=exprs(GSE129260)[, 1]),
       aes(x=expression)) +
  geom_histogram()



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------

ggplot(data=data.frame(expression_log2=log2(exprs(GSE129260)[, 1])),
       aes(x=expression_log2)) +
  geom_histogram()



## ------------------------------------------------------------------------
min(exprs(GSE129260))


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------

twoSamples <- melt(exprs(GSE129260[, 1:2]))
twoSamples$value <- log2(twoSamples$value)

ggplot(data=twoSamples, aes(x=value)) +
  facet_grid(~Var2) + geom_histogram()



## ------------------------------------------------------------------------
colSums(exprs(GSE129260))


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------

exprs(GSE129260) <- normalizeBetweenArrays(log2(exprs(GSE129260)+1), method="quantile")
twoSamples <- melt(exprs(GSE129260[, 1:2]))

ggplot(data=twoSamples, aes(x=value)) +
  facet_grid(~Var2) + geom_histogram()



## ----eval=F--------------------------------------------------------------
head(fData(GSE129260), 1000)


## ------------------------------------------------------------------------
GSE129260 <- GSE129260[!grepl("///", fData(GSE129260)$`Gene symbol`), ]
GSE129260 <- GSE129260[fData(GSE129260)$`Gene symbol` != "", ]

fData(GSE129260)$mean_expression <- apply(exprs(GSE129260), 1, mean)
GSE129260 <- GSE129260[order(fData(GSE129260)$mean_expression, decreasing = TRUE), ]
GSE129260 <- GSE129260[!duplicated(fData(GSE129260)$`Gene ID`), ]
GSE129260 <- GSE129260[seq_len(12000), ]
dim(GSE129260)


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
pcas <- prcomp(t(exprs(GSE129260)), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(GSE129260))
ggplot(plotData, aes(x=PC1, y=PC2, color=Cell, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)

mm <- exprs(GSE129260)

batch <- pData(GSE129260)$Replicate
modcombat <- model.matrix(~1, data=pData(GSE129260))
combat_mm = ComBat(dat=mm, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

pcas <- prcomp(t(combat_mm), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(GSE129260))
ggplot(plotData, aes(x=PC1, y=PC2, color=Cell, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
ggplot(plotData, aes(x=PC1, y=PC2, color=Replicate, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)



## ------------------------------------------------------------------------
fData(GSE129260)[fData(GSE129260)$`Gene symbol` == "Il10", ]


## ------------------------------------------------------------------------
exprs(GSE129260)["1450330_at", ]


## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
variance <- pcas$sdev^2
ggplot(data=data.frame(component=1:8, variance=variance),
       aes(x=component, y=variance)) +
  geom_point() + geom_line() + theme_bw()



## ----fig.height=3, fig.fullwidth=T, dev='svg', message=F-----------------
variance <- variance / sum(variance)
ggplot(data=data.frame(component=1:8, percent=variance * 100),
       aes(x=component, y=percent)) +
  geom_point() + geom_line() + theme_bw()



## ------------------------------------------------------------------------

GSE129260.design <- model.matrix(~0+Cell+Treatment+Replicate, data=pData(GSE129260))
colnames(GSE129260.design) <- c("il10neg", "il10pos", "LPS", "rep2")

fit <- lmFit(GSE129260, GSE129260.design)

fit2 <- contrasts.fit(fit, makeContrasts(il10pos - il10neg, levels=GSE129260.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## ------------------------------------------------------------------------
head(de)


## ------------------------------------------------------------------------

GSE129260.design <- model.matrix(~0+Replicate+Treatment+Cell, data=pData(GSE129260))
colnames(GSE129260.design) <- c("rep1", "rep2", "LPS", "pos")

fit <- lmFit(GSE129260, GSE129260.design)

fit2 <- contrasts.fit(fit, makeContrasts(rep2-rep1, levels=GSE129260.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## ------------------------------------------------------------------------
head(de, 100)

exprs(GSE129260)[c("1417210_at", "1426438_at"), ]

## ------------------------------------------------------------------------
## MODELING:
## ------------------------------------------------------------------------

trueAB <- matrix(c(
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
  10, 9, 8, 7, 6, 5, 4, 3, 2, 1
), ncol=2)
colnames(trueAB) <- c("A", "B")
rownames(trueAB) <- paste0("Gene ", 1:10)
head(trueAB)



## ------------------------------------------------------------------------

set.seed(1)
observed <- trueAB[, c(1, 1, 1, 2, 2, 2)]
colnames(observed) <- c("A1", "A2", "A3", "B1", "B2", "B3")
rownames(observed) <- paste0("Gene ", 1:10)
observed <- observed + rnorm(60)
head(observed)



## ------------------------------------------------------------------------

modelMatrix <- matrix(
  c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1),
  ncol = 2
)
colnames(modelMatrix) <- c("A", "B")
rownames(modelMatrix) <- colnames(observed)
head(modelMatrix)



## ------------------------------------------------------------------------

trueAB %*% t(modelMatrix)



## ------------------------------------------------------------------------

means <- observed %*% ginv(t(modelMatrix))
head(means)



## ------------------------------------------------------------------------
t.test(observed[4, 1:3], observed[4, 4:6], var.equal=TRUE)

