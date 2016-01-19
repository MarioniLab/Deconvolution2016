# This generates a plot indicating why the median fails with high numbers of low counts.

set.seed(100)
par(mfrow=c(1,2))

y1 <- rnbinom(10000, mu=20, size=0.3)
y2 <- rnbinom(10000, mu=50, size=1)

upper <- 50
pdf("small_cell.pdf")
hist(y1[y1<=upper], include.lowest=FALSE, breaks=-1:upper, main="Small cell", cex.axis=1.2, cex.main=1.4, cex.lab=1.4, xlab="Count ratio")
hist(y1[y1==0], include.lowest=FALSE, breaks=-1:upper, add=TRUE, col="black")
abline(v=median(y1), col="red")
abline(v=median(y1[y1>0]), col="red", lty=2)
dev.off()

pdf("large_cell.pdf")
hist(y2[y2<=upper], include.lowest=FALSE, breaks=-1:upper, main="Large cell", cex.axis=1.2, cex.main=1.4, cex.lab=1.4, xlab="Count ratio")
hist(y2[y2==0], include.lowest=FALSE, breaks=-1:upper, add=TRUE, col="black")
abline(v=median(y2), col="red")
abline(v=median(y2[y2>0]), col="red", lty=2)
dev.off()
