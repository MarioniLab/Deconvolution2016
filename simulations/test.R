# Defining functions.

require(AK47)

fitNorm <- function(nf, truth) {
    resids <- log2(nf) - log2(truth)
    lm(resids ~ 1, na.action=na.exclude)
}

# Generating some really noisy count data.

set.seed(100)
ngenes <- 10000
nlibs <- 200
true.facs <- 2^rnorm(nlibs, sd=0.5)
noisy.counts <- matrix(rnbinom(ngenes*nlibs, mu=10*true.facs, size=0.1), 
                       ncol=nlibs, byrow=TRUE)

# TMM with raw counts:
require(edgeR)
tmm.sf <- calcNormFactors(noisy.counts) * colSums(noisy.counts)
tmm.fitted <- fitNorm(tmm.sf, true.facs)

# TMM with averaged counts:
combined <- cbind(rowSums(noisy.counts), noisy.counts)
tmm2.sf <- calcNormFactors(combined) * colSums(combined)
tmm2.sf <- tmm2.sf[-1]
tmm2.fitted <- fitNorm(tmm2.sf, true.facs)

# Size factors with raw counts (must remove zeros for both geometric mean and for each library):
logvals <- log(noisy.counts)
logvals[is.infinite(logvals)] <- 0
gm <- exp(rowMeans(logvals))
size.sf <- apply(noisy.counts, 2, function(u) median((u/gm)[u > 0]))
size.fitted <- fitNorm(size.sf, true.facs)

# Size factors with averaged counts (must remove zeros in each library):
averaged <- rowMeans(noisy.counts)
size2.sf <- apply(noisy.counts, 2, function(u) median((u/averaged)[u > 0]))
size2.fitted <- fitNorm(size2.sf, true.facs)

# Size factors with summation:
final.sf <- normalizeBySums(noisy.counts)
final.fitted <- fitNorm(final.sf, true.facs)

pdf("benefit.pdf")
true.log.facs <- log2(true.facs)
all.range <- range(true.log.facs) + c(-0.1, 0.1)
plot(0, 0, type="n", xlim=all.range, ylim=all.range, ylab="Observed", xlab="Truth")

points(true.log.facs, true.log.facs + residuals(tmm.fitted), pch=16, col="grey")
points(true.log.facs, true.log.facs + residuals(tmm2.fitted), pch=16, col="dodgerblue")
points(true.log.facs, true.log.facs + residuals(size.fitted), pch=16, col="orangered")
points(true.log.facs, true.log.facs + residuals(size2.fitted), pch=16, col="darkred")
points(true.log.facs, true.log.facs + residuals(final.fitted), pch=16, col="black")

abline(0, 1, col="red", lwd=2)
legend("bottomright", legend=c("TMM", "TMM vs average", "Size factor (GM)", "Size factor (AM)", "Deconvolved"),
       pch=16, col=c("grey", "dodgerblue", "orangered", "darkred", "black"))
dev.off()

