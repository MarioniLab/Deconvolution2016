# This checks the brittleness of the summation approach in the presence of DE genes.

require(AK47)
require(edgeR)

set.seed(100)
ngenes <- 10000
popsize <- 250

###########################################################################

output.dir <- "results"
dir.create(output.dir, showWarning=FALSE)

make.plot <- function(sf, truth, name, main="") {
    resids <- log2(sf) - log2(truth)
    fitted <- lm(resids ~ 1, na.action=na.exclude)
    all.range <- range(truth)
    col <- rep(c(rgb(0,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0.5,0,0.5)), each=popsize)

    pdf(file.path(output.dir, paste0(name, ".pdf")))
    plot(truth, truth * 2^residuals(fitted), xlim=all.range, ylim=all.range, 
         ylab="Estimated factors", xlab="True factors", log="xy", pch=16, 
         col=col, cex.axis=1.2, cex.lab=1.4, main=main, cex.main=1.4)
    abline(0, 1, col="red")
    dev.off()
}

###########################################################################

for (scenario in 1:4) { 
    true.means <- 2^runif(ngenes, 3, 6)
    dispersions <- 2 + 100/true.means
    
    up.per.pop <- c(0.2, 0.5, 0.8)
    fc <- c(5, 5, 5)
    if (scenario==1L) { 
        nde <- 0
    } else if (scenario==2L) {
        nde <- 1000
    } else if (scenario==3L) { 
        nde <- 3000
    } else if (scenario==4L) {
        nde <- 3000
        up.per.pop <- c(0.5, 0.5, 0.5)
        fc <- c(2, 5, 10) 
    }

    counts <- list()
    true.facs <- list()
    for (x in seq_along(up.per.pop)) { 
        all.facs <- 2^rnorm(popsize, sd=0.5)
        true.facs[[x]] <- all.facs
        effective.means <- outer(true.means, all.facs, "*")
    
        chosen <- nde * (x-1) + seq_len(nde)
        is.up <- seq_len(nde*up.per.pop[x])
        upregulated <- chosen[is.up]
        downregulated <- chosen[-is.up]
        effective.means[upregulated,] <- effective.means[upregulated,] * fc[x]
        effective.means[downregulated,] <- effective.means[downregulated,] / fc[x]
    
        counts[[x]] <- matrix(rnbinom(ngenes*popsize, mu=effective.means, size=1/dispersions), ncol=popsize)
    }

    counts <- do.call(cbind, counts)
    true.facs <- unlist(true.facs)

    # TMM with raw counts:
    tmm.sf <- calcNormFactors(counts) * colSums(counts)
    make.plot(tmm.sf, true.facs, paste0("TMM_", scenario), main="TMM")

    # TMM with averaged counts:
    combined <- cbind(rowSums(counts), counts)
    tmm2.sf <- calcNormFactors(combined) * colSums(combined)
    tmm2.sf <- tmm2.sf[-1]
    make.plot(tmm2.sf, true.facs, paste0("TMMave_", scenario), main="TMM against average")

    # Size factors with raw counts (must remove zeros for both geometric mean and for each library):
    logvals <- log(counts)
    logvals[is.infinite(logvals)] <- 0
    gm <- exp(rowMeans(logvals))
    size.sf <- apply(counts, 2, function(u) median((u/gm)[u > 0]))
    make.plot(size.sf, true.facs, paste0("size_", scenario), main="Size factor")

    # Size factors with averaged counts (must remove zeros in each library):
    averaged <- rowMeans(counts)
    size2.sf <- apply(counts, 2, function(u) median((u/averaged)[u > 0]))
    make.plot(size.sf, true.facs, paste0("sizeAM_", scenario), main="Size factor with arithmetic mean")

    # Library size
    lib.sf <- colSums(counts)
    make.plot(lib.sf, true.facs, paste0("lib_", scenario), main="Library size")

    # Size factors with summation:
    final.sf <- normalizeBySums(counts, clusters=NULL)
    make.plot(final.sf, true.facs, paste0("sum_", scenario), main="Deconvolution")

    # Size factors with clustering prior to summation:
    emp.clusters <- quickCluster(counts)
    final2.sf <- normalizeBySums(counts, clusters=emp.clusters)
    make.plot(final2.sf, true.facs, paste0("sumClust_", scenario), main="Deconvolution with clustering")
}

