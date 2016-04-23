# This describes performance when pooling is random.

require(scran)
set.seed(100)

collected.order <- collected.random <- list()
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

size <- 20L
use.ave.cell <- rowMeans(exprs)
sphere <- scran:::.generateSphere(lib.sizes)
out <- scran:::.create_linear_system(ngenes, ncells, exprs, sphere, size, use.ave.cell)

design <- as.matrix(out$design)
output <- out$output
est <- solve(qr(design), output) * lib.sizes

# Trying with the opposite case, where everyone is mixed together.
sphere <- sample(ncells)
sphere <- as.integer(c(sphere, sphere))
out2 <- scran:::.create_linear_system(ngenes, ncells, exprs, sphere, size, use.ave.cell)

design2 <- as.matrix(out2$design)
output2 <- out2$output
est2 <- solve(qr(design2), output2) * lib.sizes

collected.order[[it]] <- mad(log(est/all.facs))
collected.random[[it]] <- mad(log(est2/all.facs))

cat("Ordered:", collected.order[[it]], "\n")
cat("Random:", collected.random[[it]], "\n")
cat("\n")

}

mean(unlist(collected.order))
mean(unlist(collected.random))
sd(unlist(collected.order))/sqrt(length(collected.order))
sd(unlist(collected.random))/sqrt(length(collected.random))

sessionInfo()
