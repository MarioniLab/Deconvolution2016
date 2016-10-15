# This compares the variance of the size factor estimates with and without weights. 

require(scran)
set.seed(100)

collected <- collected.w <- list()
for (it in 1:10) { 

ngenes <- 10000L
ncells <- 200L
true.means <- rgamma(ngenes, 2, 2)
dispersions <- 0.1

all.facs <- runif(ncells, 0.1, 1)
effective.means <- outer(true.means, all.facs, "*")
counts <- matrix(rnbinom(ngenes*ncells, mu=effective.means, size=1/dispersions), ncol=ncells)
lib.sizes <- colSums(counts)
exprs <- t(t(counts)/lib.sizes)

use.ave.cell <- rowMeans(exprs)
sphere <- scran:::.generateSphere(lib.sizes)
out <- scran:::.create_linear_system(exprs, sphere=sphere, sizes=c(20L, 40L, 60L, 80L, 100L), use.ave.cell)

# No weights.
design <- as.matrix(out$design)
pool.facs <- out$output
est <- solve(qr(design), pool.facs) * lib.sizes

# Precision weights based on number of cells contributing to each pool (assuming additive variances from all cells).
all.weights <- rep(1/rowSums(design))
design.w <- design * sqrt(all.weights)
pool.facs.w <- pool.facs * sqrt(all.weights)
est.w <- solve(qr(design.w), pool.facs.w) * lib.sizes

collected[[it]] <- mad(log(est/all.facs))
collected.w[[it]] <- mad(log(est.w/all.facs))

cat("No weights:", collected[[it]], "\n")
cat("Weighted:", collected.w[[it]], "\n")
cat("\n")

}

mean(unlist(collected))
mean(unlist(collected.w))
sd(unlist(collected))/sqrt(length(collected))
sd(unlist(collected.w))/sqrt(length(collected.w))

sessionInfo()
