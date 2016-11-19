# This checks the perfromance of the summation approach in the presence of increasing amounts of DE.

require(scran)
require(edgeR)
require(DESeq2)

set.seed(100)
ngenes <- 10000
popsize <- 200

###########################################################################

output.dir <- "results_de"
dir.create(output.dir, showWarning=FALSE)

make.plot <- function(sf, truth, main="") {
    resids <- log2(sf) - log2(truth)
    fitted <- lm(resids ~ 1, na.action=na.exclude)
    all.range <- range(truth)
    plot(truth, truth * 2^residuals(fitted), xlim=all.range, ylim=all.range, 
         ylab="Estimated factors", xlab="True factors", log="xy", pch=16, 
         cex.axis=1.5, cex.lab=1.8, main=main, cex.main=1.8)
    abline(0, 1, col="red")
}

###########################################################################

pdf(file.path(output.dir, "comparator.pdf"), width=10, height=6)
par(mfrow=c(1,2))
for (num in c(10L, 50L, 100L, 500L, 1000L, 2000L)){ 
    true.means <- rgamma(ngenes, 2, 2)
    dispersions <- 0.1
    
    all.facs <- 2^rnorm(popsize, sd=0.5)
    effective.means <- outer(true.means, all.facs, "*")
    effective.means[1:num,1:100] <- effective.means[1:num,1:100]*10 # Upregulation of a subset of DE genes.

    counts <- matrix(rnbinom(ngenes*popsize, mu=effective.means, size=1/dispersions), ncol=popsize)
    true.facs <- all.facs

    # Size factors with summation:
    final.sf <- computeSumFactors(counts, clusters=NULL, sizes=c(20, 40, 60, 80, 100))
    make.plot(final.sf, true.facs, main=sprintf("Deconvolution (%.1f%% DE)", num/ngenes*100))

    lib.sf <- colSums(counts)
    make.plot(lib.sf, true.facs, main=sprintf("Library size (%.1f%% DE)", num/ngenes*100))
}
dev.off()
