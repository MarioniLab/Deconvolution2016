# This checks the brittleness of the summation approach in the presence of DE genes.

require(AK47)
require(edgeR)

set.seed(100)
ngenes <- 10000
popsize <- 200
subpops <- 3
nlibs <- popsize * subpops

###########################################################################

output.dir <- "results"
dir.create(output.dir, showWarning=FALSE)

make.plot <- function(sf, truth, name) {
    resids <- log2(sf) - log2(truth)
    fitted <- lm(resids ~ 1, na.action=na.exclude)
    all.range <- range(truth)

    pdf(file.path(output.dir, paste0(name, ".pdf")))
    plot(truth, truth * 2^residuals(fitted), xlim=all.range, ylim=all.range, 
         ylab="Estimated factors", xlab="True factors", log="xy", pch=16, 
         col=rgb(0,0,0,0.5), cex.axis=1.2, cex.lab=1.4)
    abline(0, 1, col="red")
    dev.off()
}

###########################################################################

for (scenario in 1:3) { 
    true.facs <- 2^rnorm(nlibs, sd=0.5)
    true.means <- 2^runif(ngenes, 3, 6)
    effective.means <- outer(true.means, true.facs, "*")
    dispersions <- 2 + 100/true.means
    
    if (scenario==1L) { 
        clusters <- NULL
    } else {
        # Each subpopulation has its own DE set, which is upregulated to a different extent.
        if (scenario==2L) {
            nde <- 1000
            clusters <- rep(seq_len(subpops), each=popsize)
        } else if (scenario==3L) { 
            nde <- 3000
            clusters <- rep(seq_len(subpops), each=popsize)
        }

        last.gene <- last.cell <- 0L
        for (x in seq_len(subpops)) { 
            curgene <- last.gene + seq_len(nde)
            curcell <- last.cell + seq_len(popsize)
            effective.means[curgene,curcell] <- effective.means[curgene,curcell] * (x+2) # unique FC per group.
            last.gene <- last.gene + nde
            last.cell <- last.cell + popsize
        }
    }
    counts <- matrix(rnbinom(ngenes*nlibs, mu=effective.means, size=1/dispersions), ncol=nlibs)

    # TMM with raw counts:
    tmm.sf <- calcNormFactors(counts) * colSums(counts)
    make.plot(tmm.sf, true.facs, paste0("TMM_", scenario))

    # TMM with averaged counts:
    combined <- cbind(rowSums(counts), counts)
    tmm2.sf <- calcNormFactors(combined) * colSums(combined)
    tmm2.sf <- tmm2.sf[-1]
    make.plot(tmm2.sf, true.facs, paste0("TMMave_", scenario))

    # Size factors with raw counts (must remove zeros for both geometric mean and for each library):
    logvals <- log(counts)
    logvals[is.infinite(logvals)] <- 0
    gm <- exp(rowMeans(logvals))
    size.sf <- apply(counts, 2, function(u) median((u/gm)[u > 0]))
    make.plot(size.sf, true.facs, paste0("size_", scenario))

    # Size factors with averaged counts (must remove zeros in each library):
    averaged <- rowMeans(counts)
    size2.sf <- apply(counts, 2, function(u) median((u/averaged)[u > 0]))
    make.plot(size.sf, true.facs, paste0("sizeAM_", scenario))

    # Size factors with summation:
    final.sf <- normalizeBySums(counts, clusters=NULL)
    make.plot(final.sf, true.facs, paste0("sum_", scenario))

    # Size factors with clustering prior to summation:
    if (scenario > 1L) {
        final2.sf <- normalizeBySums(counts, clusters=clusters)
        make.plot(final2.sf, true.facs, paste0("sumClust_", scenario))
    }
}

