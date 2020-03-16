
## Installing packages if they are not installed

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("affxparser", quietly = TRUE)) BiocManager::install("affxparser")
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("mouse4302.db", quietly = TRUE)) BiocManager::install("mouse4302.db")



## Reading one cell file (as an example)

library(affxparser)
library(affy)

CELfile <- readCel("GSE129260_RAW/GSM3703675_IL-10_posi_anti-CD40-1.CEL")



## ------------------------------------------------------------------------
head(CELfile$header)


## ------------------------------------------------------------------------
head(CELfile$intensities)


## Reading all CEL files and performing RMA normalization

files <- list.files("GSE129260_RAW/", full.names = T)
microarrayData <- justRMA(filenames = files)


## Showing some probe-level expression values
exprs(microarrayData)[1:5, 1:2]


## Symbol/probe annotation
library(mouse4302.db)

symbolAnnotation <- as.list(mouse4302SYMBOL)
head(symbolAnnotation, 3)



## Getting the dataset from GEO

library(GEOquery)
GSE129260 <- getGEO("GSE129260", AnnotGPL = TRUE)[[1]]



## Dataset dimensions
dim(exprs(GSE129260))



## Some dataset values

head(exprs(GSE129260))



## Phenotypical data

head(pData(GSE129260)[, 1:2])



## Feature data

head(fData(GSE129260)[, 1:2])


