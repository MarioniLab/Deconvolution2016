# This calculates the standard error of the size factor estimates from each method.

require(scran)
require(edgeR)
require(DESeq2)

set.seed(100)
ngenes <- 10000
popsize <- 250

###########################################################################

recollect <- function(sf, collected) {
    sf <- log(sf)
    collected[[length(collected)+1]] <- exp(sf - mean(sf))    
    return(collected)
}

###########################################################################

true.means <- rgamma(ngenes, 2, 2)
dispersions <- 0.1

up.per.pop <- c(0.2, 0.5, 0.8)
fc.up <- c(5, 5, 5)
fc.down <- c(0, 0, 0)
nde <- 0 

collected.tmm <- collected.deseq <- collected.lib <- list()
collected.deconv <- collected.deconv.more <- collected.deconv.extra <- collected.large <- collected.small <- list()      
all.facs <- 2^rnorm(popsize, sd=0.5) # Constant factors!

for (it in 1:10) {
    counts <- list()
    true.facs <- list()
    for (x in seq_along(up.per.pop)) { 
        true.facs[[x]] <- all.facs
        effective.means <- outer(true.means, all.facs, "*")
    
        chosen <- nde * (x-1) + seq_len(nde)
        is.up <- seq_len(nde*up.per.pop[x])
        upregulated <- chosen[is.up]
        downregulated <- chosen[-is.up]
        effective.means[upregulated,] <- effective.means[upregulated,] * fc.up[x]
        effective.means[downregulated,] <- effective.means[downregulated,] * fc.down[x]
    
        counts[[x]] <- matrix(rnbinom(ngenes*popsize, mu=effective.means, size=1/dispersions), ncol=popsize)
    }

    counts <- do.call(cbind, counts)
    true.facs <- unlist(true.facs)

    # TMM with raw counts:
    tmm.sf <- calcNormFactors(counts) * colSums(counts)
    collected.tmm <- recollect(tmm.sf, collected.tmm)

    # Size factors with raw counts (must counter zeros for both geometric mean and for each library):
    logvals <- log(counts)
    logvals[is.infinite(logvals)] <- NA_real_
    gm <- exp(rowMeans(logvals, na.rm=TRUE))
    size.sf <- estimateSizeFactorsForMatrix(counts, geoMeans=gm)
    collected.deseq <- recollect(size.sf, collected.deseq)
   
    # Library size
    lib.sf <- colSums(counts)
    collected.lib <- recollect(lib.sf, collected.lib)

    # Size factors with clustering prior to summation:
    emp.clusters <- quickCluster(counts)
    final2.sf <- computeSumFactors(counts, clusters=emp.clusters, sizes=1:5*20)
    collected.deconv <- recollect(final2.sf, collected.deconv)
    final3.sf <- computeSumFactors(counts, clusters=emp.clusters, sizes=2:10*10)
    collected.deconv.more <- recollect(final3.sf, collected.deconv.more)
    final4.sf <- computeSumFactors(counts, clusters=emp.clusters, sizes=10:50*2)
    collected.deconv.extra <- recollect(final4.sf, collected.deconv.extra)

    # Size factors using only small or large pool sizes.
    final.large.sf <- suppressWarnings(computeSumFactors(counts, size=200, clusters=emp.clusters))
    collected.large <- recollect(final.large.sf, collected.large)
    final.small.sf <- suppressWarnings(computeSumFactors(counts, size=20, clusters=emp.clusters))
    collected.small <- recollect(final.small.sf, collected.small)
}

# Precision of estimates.
summary(apply(log(do.call(rbind, collected.tmm)), 2, mad))
summary(apply(log(do.call(rbind, collected.deseq)), 2, mad))
summary(apply(log(do.call(rbind, collected.lib)), 2, mad))

summary(apply(log(do.call(rbind, collected.deconv)), 2, mad))
summary(apply(log(do.call(rbind, collected.deconv.more)), 2, mad))
summary(apply(log(do.call(rbind, collected.deconv.extra)), 2, mad))

summary(apply(log(do.call(rbind, collected.large)), 2, mad))
summary(apply(log(do.call(rbind, collected.small)), 2, mad))
mad(log(all.facs))
