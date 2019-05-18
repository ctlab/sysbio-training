library(tximport)

library(data.table)
tx2gene <- fread("/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.tx2gene.tsv")

kallisto.files <- list.files("./kallisto", pattern="abundance.tsv", recursive = TRUE, full.names = TRUE)

txi.kallisto.tsv <- tximport(files=kallisto.files, 
                             type = "kallisto", 
                             tx2gene = tx2gene, ignoreAfterBar = TRUE)
head(txi.kallisto.tsv$counts)

fc.files <- list.files("./featureCounts/", pattern="fc.txt$", recursive = TRUE, full.names = TRUE)
fc_res <- lapply(fc.files, fread)

fc_mat <- do.call(cbind, lapply(fc_res, function(x) x[[ncol(x)]]))
rownames(fc_mat) <- fc_res[[1]][, Geneid]

nrow(fc_mat)
nrow(txi.kallisto.tsv$counts)

kallisto_v <- setNames(txi.kallisto.tsv$counts[,1], rownames(txi.kallisto.tsv$counts))
fc_v <- setNames(fc_mat[,1], rownames(fc_mat))
kallisto_v <- kallisto_v[names(fc_v)]
plot(log2(fc_v+1), log2(kallisto_v+1))
