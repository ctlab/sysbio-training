#!/usr/bin/env Rscript
library(GenomicFeatures)
library(data.table)
TxDb <- makeTxDbFromGFF(file = "gencode.v29.annotation.gtf")

tx2gene <- select(TxDb, keys(TxDb, keytype = "TXNAME"), column="GENEID", keytype="TXNAME")
fwrite(tx2gene, file="gencode.v29.annotation.tx2gene.tsv", sep="\t")

