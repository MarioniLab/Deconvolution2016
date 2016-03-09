# This describes performance when pooling is random.

require(scran)
set.seed(100)

for (it in 1:10) { 

ngenes <- 10000L
ncells <- 200L
true.means <- 2^runif(ngenes, 0, 3)
dispersions <- 1

all.facs <- runif(ncells, 0.1, 1)
effective.means <- outer(true.means, all.facs, "*")
counts <- matrix(rnbinom(ngenes*ncells, mu=effective.means, size=1/dispersions), ncol=ncells)
lib.sizes <- colSums(counts)
exprs <- t(t(counts)/lib.sizes)

size <- 20L
use.ave.cell <- rowMeans(exprs)
sphere <- scran:::.generateSphere(lib.sizes) - 1L
out <- .Call(scran:::cxx_forge_system, ngenes, ncells, exprs, sphere, size, use.ave.cell)

out.ex <- .Call(scran:::cxx_forge_system, ngenes, ncells, exprs, sphere, 1L, use.ave.cell)
design <- c(out[1], out.ex[1])
output <- c(out[2], out.ex[2])
design <- do.call(rbind, design)
output <- unlist(output)

weights <- rep(c(1, 0.00001), c(nrow(design)-ncells, ncells))
root.weights <- sqrt(weights)
design <- design * root.weights
output <- output * root.weights
est <- solve(qr(design), output) * lib.sizes

# Trying with the opposite case, where everyone is mixed together.
sphere <- sample(ncells)
sphere <- as.integer(c(sphere, sphere) - 1)
out2 <- .Call(scran:::cxx_forge_system, ngenes, ncells, exprs, sphere, size, use.ave.cell)
design2 <- c(out2[1], out.ex[1])
output2 <- c(out2[2], out.ex[2])
design2 <- do.call(rbind, design2)
output2 <- unlist(output2)
design2 <- design2 * root.weights
output2 <- output2 * root.weights
est2 <- solve(qr(design2), output2) * lib.sizes

cat("Ordered:", mad(log(est/all.facs)), "\n")
cat("Random:", mad(log(est2/all.facs)), "\n")
cat("\n")

}

sessionInfo()
