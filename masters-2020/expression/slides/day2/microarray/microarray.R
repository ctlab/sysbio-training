## ------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("affxparser", quietly = TRUE)) BiocManager::install("affxparser")
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("mouse4302.db", quietly = TRUE)) BiocManager::install("mouse4302.db")



## ----cache=TRUE----------------------------------------------------------
library(affxparser)
library(affy)

CELfile <- readCel("GSE129260_RAW/GSM3703675_IL-10_posi_anti-CD40-1.CEL")



## ------------------------------------------------------------------------
head(CELfile$header)


## ------------------------------------------------------------------------
head(CELfile$intensities)


## ----cache=TRUE, suppressMessages=TRUE-----------------------------------
files <- list.files("GSE129260_RAW/", full.names = T)
microarrayData <- justRMA(filenames = files)


## ----cache=TRUE----------------------------------------------------------
exprs(microarrayData)[1:5, 1:2]


## ----message=F-----------------------------------------------------------
library(mouse4302.db)

symbolAnnotation <- as.list(mouse4302SYMBOL)
head(symbolAnnotation, 3)



## ----message=F, cache=T--------------------------------------------------

library(GEOquery)
GSE129260 <- getGEO("GSE129260", AnnotGPL = TRUE)[[1]]



## ----message=F-----------------------------------------------------------
dim(exprs(GSE129260))
head(exprs(GSE129260))



## ----message=F-----------------------------------------------------------

head(pData(GSE129260)[, 1:2])



## ----message=F-----------------------------------------------------------

head(fData(GSE129260)[, 1:2])


