## Installing packages (in case you don't have them)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

library(GEOquery)
library(Biobase)
library(ggplot2)
library(reshape2)
library(limma)
library(MASS)


## Let's load the dataset

GSE129260 <- getGEO("GSE129260", AnnotGPL = TRUE)[[1]]



## Keeping only relevant phenotypical data
pData(GSE129260)$rep <- gsub(".*(rep\\d)$", "\\1", pData(GSE129260)$title)
pData(GSE129260) <- pData(GSE129260)[, c("characteristics_ch1.1", "characteristics_ch1.2", "rep")]

colnames(pData(GSE129260)) <- c("Cell", "Treatment", "Replicate")
head(pData(GSE129260))



## Feature data
colnames(fData(GSE129260))



## Keeping only relevant feature data
fData(GSE129260) <- fData(GSE129260)[, c("ID", "Gene symbol", "Gene ID")]
head(fData(GSE129260))



## Looking at the original distribution

ggplot(data=data.frame(expression=exprs(GSE129260)[, 1]),
       aes(x=expression)) +
  geom_histogram()



## Looking at log-distribution

ggplot(data=data.frame(expression_log2=log2(exprs(GSE129260)[, 1])),
       aes(x=expression_log2)) +
  geom_histogram()



## minimal/maximal value is also an indication
min(exprs(GSE129260))
max(exprs(GSE129260))


## Example of log-tranformed (two samples next to each other)

twoSamples <- melt(exprs(GSE129260[, 1:2]))
twoSamples$value <- log2(twoSamples$value)

ggplot(data=twoSamples, aes(x=value)) +
  facet_grid(~Var2) + geom_histogram()



## Actual normalization
## normalizeBetweenArrays - quantile normalization (method="quantile")
## please note (log2(data + 1)) - we are log-transforming data before quantile normalizing it

exprs(GSE129260) <- normalizeBetweenArrays(log2(exprs(GSE129260)+1), method="quantile")

## Two samples after transformation

twoSamples <- melt(exprs(GSE129260[, 1:2]))

ggplot(data=twoSamples, aes(x=value)) +
  facet_grid(~Var2) + geom_histogram()


## Mapping probes to genes:
## 1. removing probes mapped to many genes
## 2. removing probes mapper to no genes
## 3. calculating average probe expression
## 4. sorting all probes by average expression
## 5. For each gene we only keep most-expressed probes
## 6. leaving only top 12 000 expressed genes

GSE129260 <- GSE129260[!grepl("///", fData(GSE129260)$`Gene symbol`), ]
GSE129260 <- GSE129260[fData(GSE129260)$`Gene symbol` != "", ]

fData(GSE129260)$mean_expression <- apply(exprs(GSE129260), 1, mean)
GSE129260 <- GSE129260[order(fData(GSE129260)$mean_expression, decreasing = TRUE), ]
GSE129260 <- GSE129260[!duplicated(fData(GSE129260)$`Gene ID`), ]
GSE129260 <- GSE129260[seq_len(12000), ]
dim(GSE129260)


## Showing PCA
pcas <- prcomp(t(exprs(GSE129260)), scale. = T)
plotData <- cbind(pcas$x[, 1:2], pData(GSE129260))
ggplot(plotData, aes(x=PC1, y=PC2, color=Cell, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)



## Showing PCA with different kinds of legend
ggplot(plotData, aes(x=PC1, y=PC2, color=Replicate, shape=Treatment)) +
  geom_point() + theme_bw() + theme(aspect.ratio = 1)



## Finding probe with Il10
fData(GSE129260)[fData(GSE129260)$`Gene symbol` == "Il10", ]


## Checking Il10 expression
exprs(GSE129260)["1450330_at", ]


## Showing variance (calculated with PCA)
variance <- pcas$sdev^2
ggplot(data=data.frame(component=1:8, variance=variance),
       aes(x=component, y=variance)) +
  geom_point() + geom_line() + theme_bw()



## Showing % of variance (calculated with PCA)
variance <- variance / sum(variance)
ggplot(data=data.frame(component=1:8, percent=variance * 100),
       aes(x=component, y=percent)) +
  geom_point() + geom_line() + theme_bw()



## Differential expression (Il10pos vs Il10neg)

GSE129260.design <- model.matrix(~0+Cell+Treatment+Replicate, data=pData(GSE129260))
colnames(GSE129260.design) <- c("il10neg", "il10pos", "LPS", "rep2")

fit <- lmFit(GSE129260, GSE129260.design)

fit2 <- contrasts.fit(fit, makeContrasts(il10pos - il10neg, levels=GSE129260.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## ------------------------------------------------------------------------
head(de)


## Differential expression (Replicate 1 vs Replicate 2)

GSE129260.design <- model.matrix(~0+Replicate+Treatment+Cell, data=pData(GSE129260))
colnames(GSE129260.design) <- c("rep1", "rep2", "LPS", "pos")

fit <- lmFit(GSE129260, GSE129260.design)

fit2 <- contrasts.fit(fit, makeContrasts(rep2-rep1, levels=GSE129260.design))
fit2 <- eBayes(fit2, trend = T)

de <- topTable(fit2, adjust.method="BH", number=Inf, sort.by = "P")



## ------------------------------------------------------------------------
head(de)

