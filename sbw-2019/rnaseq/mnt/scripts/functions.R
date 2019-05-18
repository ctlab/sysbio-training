#' Makes nice PCA plot for expression data
#' @param es an ExpressionSet object, should be normalized
#' @param c1 a number of the first component to plot (numeric)
#' @param c2 a number of the second component to plot (numeric)
#' @examples 
#' pcaPlot(es.norm, 1, 2) + aes(color=time)
# from https://github.com/assaron/r-utils/blob/master/R/exprs.R
pcaPlot <- function(es, c1, c2) {
  stopifnot(require(ggplot2))
  pca <- prcomp(t(exprs(es)))
  
  explained <- (pca$sdev)^2 / sum(pca$sdev^2)
  
  xs <- sprintf("PC%s", seq_along(explained))
  xlabs <- sprintf("%s (%.1f%%)", xs, explained * 100)
  
  d <- cbind(as.data.frame(pca$x), pData(es))
  if (!"sample" %in% colnames(d)) {
    d <- cbind(d, sample=colnames(es))
  }
  
  pp <- ggplot(data=d)
  
  pp + aes_string(x=xs[c1], y=xs[c2]) +
    geom_point(size=3) +
    xlab(xlabs[c1]) + ylab(xlabs[c2])
}

# from https://github.com/assaron/r-utils/blob/master/R/exprs.R
write.gct <- function(es, file, gzip=FALSE) {
  stopifnot(require(Biobase))
  if (gzip) {
    con <- gzfile(file)
  } else {
    con <- file(file)
  }
  open(con, open="w")
  writeLines("#1.3", con)
  ann.col <- ncol(pData(es))
  ann.row <- ncol(fData(es))
  writeLines(sprintf("%s\t%s\t%s\t%s", nrow(es), ncol(es), ann.row, ann.col), con)
  writeLines(paste0(c("ID", colnames(fData(es)), colnames(es)), collapse="\t"), con)
  
  ann.col.table <- t(as.matrix(pData(es)))    
  ann.col.table <- cbind(matrix(rep(NA, ann.row*ann.col), nrow=ann.col), ann.col.table)
  write.table(ann.col.table, file=con, quote=F, sep="\t", row.names=T, col.names=F)                          
  write.table(cbind(fData(es), exprs(es)), file=con, quote=F, sep="\t", row.names=T, col.names=F)                          
  close(con)    
}

# from https://github.com/assaron/r-utils/blob/master/R/exprs.R
collapseBy <- function(es, factor, FUN=median) {
  ranks <- apply(exprs(es), 1, FUN)
  t <- data.frame(f=factor, i=seq_along(ranks), r=ranks)
  t <- t[order(t$r, decreasing=T), ]
  keep <- t[!duplicated(t$f) & !is.na(t$f),]$i
  res <- es[keep, ]
  fData(res)$origin <- rownames(res)
  rownames(res) <- factor[keep]
  res
}

# from https://github.com/assaron/r-utils/blob/master/R/exprs.R
#' Make contrast DESeq2-style
#' @param contrast vector of three elements (factor, level2, level1)
#' @return  contrast table for level2-level1 comparison
#' @examples 
#' es.design <- model.matrix(~0+condition, data=pData(es.norm))
#' fit <- lmFit(es.norm, es.design)
#' fit2 <- contrasts.fit(fit, makeContrasts2(c("condition", "Treatment", "Control"), 
#'                                          levels=es.design))
#' fit2 <- eBayes(fit2)
#' de <- topTable(fit2, adjust.method="BH", number=Inf)
#' @export
makeContrasts2 <- function(contrast, levels) {
  f <- contrast[1]
  c1 <- contrast[3]
  c2 <- contrast[2]
  tag <- sprintf("%s%s-%s%s", f, c2, f, c1)
  ls <- colnames(levels)
  res <- matrix(rep(0, length(ls)), nrow=length(ls), dimnames = list(Levels=ls, Contrasts=tag))
  res[paste0(f, c2), ] <- 1
  res[paste0(f, c1), ] <- -1
  res
}