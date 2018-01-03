# This checks the brittleness of the summation approach in the presence of DE genes.

require(scran)
require(edgeR)
require(DESeq2)

set.seed(100)
ngenes <- 10000
popsize <- 250

###########################################################################

output.dir <- "results_high"
dir.create(output.dir, showWarning=FALSE)

make.plot <- function(sf, truth, name, main="") {
    resids <- log2(sf) - log2(truth)
    fitted <- lm(resids ~ 1, na.action=na.exclude)
    all.range <- range(truth)
    col <- rep(c("black", "dodgerblue", "orange"), each=popsize)
    shuffle <- as.vector(t(matrix(seq_along(sf), nrow=popsize)))

    pdf(file.path(output.dir, paste0(name, ".pdf")))
    par(mar=c(5.1,5.1,4.1,1.1))
    plot(truth[shuffle], (truth * 2^residuals(fitted))[shuffle], xlim=all.range, ylim=all.range, 
         ylab="Estimated factors", xlab="True factors", log="xy", pch=16, 
         col=col[shuffle], cex.axis=1.5, cex.lab=1.8, main=main, cex.main=1.8)
    abline(0, 1, col="red")
    dev.off()
}

###########################################################################

for (scenario in 1:4) { 
    true.means <- 2^runif(ngenes, 9, 12)
    dispersions <- 2 + 100/true.means
    
    up.per.pop <- c(0.2, 0.5, 0.8)
    fc.up <- c(5, 5, 5)
    fc.down <- c(0, 0, 0)
    if (scenario==1L) { # No DE
        nde <- 0 
    } else if (scenario==2L) { # Some DE; same number of genes, different numbers per direction
        nde <- 1000
    } else if (scenario==3L) { # More DE; same number of genes, different numbers per direction
        nde <- 3000
    } else if (scenario==4L) { # More DE; same numbers per direction, different magnitude of DE
        nde <- 3000
        up.per.pop <- c(0.5, 0.5, 0.5)
        fc.up <- c(2, 5, 10) 
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
        effective.means[upregulated,] <- effective.means[upregulated,] * fc.up[x]
        effective.means[downregulated,] <- effective.means[downregulated,] * fc.down[x]
    
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

    # Size factors with raw counts (must counter zeros for both geometric mean and for each library):
    logvals <- log(counts)
    logvals[is.infinite(logvals)] <- NA_real_
    gm <- exp(rowMeans(logvals, na.rm=TRUE))
    size.sf <- estimateSizeFactorsForMatrix(counts, geoMeans=gm)
    make.plot(size.sf, true.facs, paste0("size_", scenario), main="DESeq")

    # Size factors with averaged counts (still removes zeros in each library):
    size2.sf <- estimateSizeFactorsForMatrix(counts, geoMeans=rowMeans(counts))
    make.plot(size2.sf, true.facs, paste0("sizeAM_", scenario), main="DESeq with arithmetic mean")

    # Size factors with an added pseudo-count.
    sizeP.sf <- estimateSizeFactorsForMatrix(counts+1)
    make.plot(sizeP.sf, true.facs, paste0("sizeP_", scenario), main="DESeq with pseudo-count")
    
    # Size factors with a library size-adjusted pseudo-count.
    lib.size <- colSums(counts)
    pcounts <- t(t(counts) + lib.size/mean(lib.size))
    sizeP2.sf <- estimateSizeFactorsForMatrix(pcounts)
    make.plot(sizeP2.sf, true.facs, paste0("sizeP2_", scenario), main="DESeq with library size-adjusted pseudo-count")
    
    # Library size
    lib.sf <- colSums(counts)
    make.plot(lib.sf, true.facs, paste0("lib_", scenario), main="Library size")

    # Size factors with summation:
    final.sf <- computeSumFactors(counts, clusters=NULL, min.mean=0)
    make.plot(final.sf, true.facs, paste0("sum_", scenario), main="Deconvolution")

    # Size factors with clustering prior to summation:
    emp.clusters <- quickCluster(counts)
    final2.sf <- computeSumFactors(counts, clusters=emp.clusters, min.mean=0)
    make.plot(final2.sf, true.facs, paste0("sumClust_", scenario), main="Deconvolution with clustering")
}

