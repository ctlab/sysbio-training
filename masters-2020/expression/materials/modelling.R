## Expression in true A vs B

trueAB <- matrix(c(
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
  10, 9, 8, 7, 6, 5, 4, 3, 2, 1
), ncol=2)
colnames(trueAB) <- c("A", "B")
rownames(trueAB) <- paste0("Gene ", 1:10)
head(trueAB)



## Making 3 replicates of both A and B + adding noise

set.seed(1)
observed <- trueAB[, c(1, 1, 1, 2, 2, 2)]
colnames(observed) <- c("A1", "A2", "A3", "B1", "B2", "B3")
rownames(observed) <- paste0("Gene ", 1:10)
observed <- observed + rnorm(60)
head(observed)



## Model matrix describes samples

modelMatrix <- matrix(
  c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1),
  ncol = 2
)
colnames(modelMatrix) <- c("A", "B")
rownames(modelMatrix) <- colnames(observed)
head(modelMatrix)



## Matrix multiplication will give us something close to observed matrix (but without noise)

trueAB %*% t(modelMatrix)



## We could get means in groups using multiplication with inverse matrix

means <- observed %*% ginv(t(modelMatrix))
head(means)



## T-test for one gene
t.test(observed[4, 1:3], observed[4, 4:6], var.equal=TRUE)
