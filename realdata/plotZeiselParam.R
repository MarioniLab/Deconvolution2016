cur.data <- readRDS("KleinData.rds")
counts <- cur.data$Counts
size.facs <- cur.data$SF
grouping <- cur.data$Cells

pdf("Klein_Means.pdf")
par(mar=c(5.1, 5.1, 1.1, 1.1))
rwm <- rowMeans(counts)
hist(rwm[rwm < 10], freq=FALSE, xlab="Average count per gene", col="grey80", cex.lab=1.4, cex.axis=1.2, main="")
lines(density(rgamma(1e5, 2, 2)), col="red", lwd=2)
dev.off()

require(edgeR)
y <- DGEList(counts)
y$samples$norm.factors <- size.facs$Deconvolution/y$samples$lib.size
design <- model.matrix(~grouping$Time)
y <- estimateDisp(y, design, prior.df=0, trend='none')

pdf("Klein_Dispersions.pdf")
par(mar=c(5.1, 5.1, 1.1, 1.1))
smoothScatter(log2(rwm), log2(y$tagwise.dispersion), xlab=expression(Log[2]~"average count per gene"), ylab=expression(Log[2]~"NB dispersion estimate"), cex.axis=1.2, cex.lab=1.4)
o <- order(rwm)
lines(log2(rwm)[o], log2(y$trended.dispersion)[o], col="red", lwd=2, lty=2)
abline(h=log2(0.1), col="red", lwd=2)
dev.off()

