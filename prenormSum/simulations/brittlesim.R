# This checks the brittleness of the summation approach in the presence of DE genes.

set.seed(100)
ngenes <- 10000
nlibs <- 800
require(AK47)
require(edgeR)

fitNorm <- function(nf, truth) {
    resids <- log2(nf) - log2(truth)
    lm(resids ~ 1, na.action=na.exclude)
}
            
popsize <- 200
subpops <- 3

###########################################################################

pdf("results.pdf")

for (scenario in 1:5) { 
    true.facs <- 2^rnorm(nlibs, sd=0.5)
    true.means <- 2^runif(ngenes, 3, 6)
    effective.means <- outer(true.means, true.facs, "*")
    dispersions <- 2 + 100/true.means
    
    if (scenario==1L) { 
        main <- "No DE"
        clusters <- NULL
    } else {
        # Each subpopulation has its own DE set, which is upregulated to a different extent.
        if (scenario==2L) {
            nde <- 1000
            clusters <- NULL
            extra <- ""
        } else if (scenario==3L) {
            nde <- 3000
            clusters <- NULL
            extra <- ""
        } else if (scenario==4L) { 
            nde <- 1000
            clusters <- rep(1:(subpops+1L), each=popsize)
            extra <- " (by cluster)"
        } else if (scenario==5L) { 
            nde <- 3000
            clusters <- rep(1:(subpops+1L), each=popsize)
            extra <- " (by cluster)"
        }

        last.gene <- last.cell <- 0L
        for (x in seq_len(subpops)) { 
            curgene <- last.gene + seq_len(nde)
            curcell <- last.cell + seq_len(popsize)
            effective.means[curgene,curcell] <- effective.means[curgene,curcell] * (x+2) # unique FC per group.
            last.gene <- last.gene + nde
            last.cell <- last.cell + popsize
        }
        main <- sprintf("DE in %i genes for %i subpopulations%s", nde, subpops, extra)
    }
    counts <- matrix(rnbinom(ngenes*nlibs, mu=effective.means, size=1/dispersions), ncol=nlibs)

    # TMM with raw counts:
    tmm.sf <- calcNormFactors(counts) * colSums(counts)
    tmm.fitted <- fitNorm(tmm.sf, true.facs)

    # TMM with averaged counts:
    combined <- cbind(rowSums(counts), counts)
    tmm2.sf <- calcNormFactors(combined) * colSums(combined)
    tmm2.sf <- tmm2.sf[-1]
    tmm2.fitted <- fitNorm(tmm2.sf, true.facs)

    # Size factors with raw counts (must remove zeros for both geometric mean and for each library):
    logvals <- log(counts)
    logvals[is.infinite(logvals)] <- 0
    gm <- exp(rowMeans(logvals))
    size.sf <- apply(counts, 2, function(u) median((u/gm)[u > 0]))
    size.fitted <- fitNorm(size.sf, true.facs)

    # Size factors with averaged counts (must remove zeros in each library):
    averaged <- rowMeans(counts)
    size2.sf <- apply(counts, 2, function(u) median((u/averaged)[u > 0]))
    size2.fitted <- fitNorm(size2.sf, true.facs)

    # Size factors with summation:
    final.sf <- normalizeBySums(counts, clusters=clusters)
    final.fitted <- fitNorm(final.sf, true.facs)

    # Making a plot of the results.
    true.log.facs <- log2(true.facs)
    all.range <- range(true.log.facs) + c(-0.1, 0.1)
    plot(0, 0, type="n", xlim=all.range, ylim=all.range, ylab="Observed", xlab="Truth", main=main)

    points(true.log.facs, true.log.facs + residuals(tmm.fitted), pch=16, col="grey")
    points(true.log.facs, true.log.facs + residuals(tmm2.fitted), pch=16, col="dodgerblue")
    points(true.log.facs, true.log.facs + residuals(size.fitted), pch=16, col="orangered")
    points(true.log.facs, true.log.facs + residuals(size2.fitted), pch=16, col="darkred")
    points(true.log.facs, true.log.facs + residuals(final.fitted), pch=16, col="black")

    abline(0, 1, col="red", lwd=2)
    legend("bottomright", legend=c("TMM", "TMM vs average", "Size factor (GM)", "Size factor (AM)", "Deconvolved"),
                  pch=16, col=c("grey", "dodgerblue", "orangered", "darkred", "black"))
}

dev.off()
