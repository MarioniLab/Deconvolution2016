# This checks the time complexity of the computeSumFactors.

require(scran)
set.seed(100)        
ngenes <- 10000
all.sizes <- 200*(1:4)
all.times.raw <- all.times.clust <- list()

for (it in 1:10) {
    cur.time.raw <- cur.time.clust <- NULL
    for (s in all.sizes) {
        counts <- matrix(rpois(ngenes*s, lambda=10), ncol=s)
        times <- system.time(out <- computeSumFactors(counts))
        cur.time.raw <- c(cur.time.raw, times[1])

        n <- s/200
        clusters <- rep(seq_len(n), each=200)
        times2 <- system.time(out <- computeSumFactors(counts, cluster=clusters))
        cur.time.clust <- c(cur.time.clust, times2[1])
    }

    all.times.raw[[it]] <- cur.time.raw
    all.times.clust[[it]] <- cur.time.clust
}

pdf("timings.pdf")
par(mar=c(5.1, 5.1, 1.1, 1.1))
plot(0,0, type="n", xlab="Number of cells", ylab="Computation time (seconds)", xlim=range(all.sizes), ylim=c(0, 10), cex.axis=1.2, cex.lab=1.4)
for (it in seq_along(all.times.raw)) {
    points(all.sizes, all.times.raw[[it]], pch=16, col=rgb(1,0,0,0.2))
    lines(all.sizes, all.times.raw[[it]], col=rgb(1,0,0,0.2), lwd=2)
    points(all.sizes, all.times.clust[[it]], pch=16, col=rgb(0,0,1,0.2))
    lines(all.sizes, all.times.clust[[it]], col=rgb(0,0,1,0.2), lwd=2)
}
legend("bottomright", lwd=2, col=c("red", "blue"), legend=c("Without clustering", "With clustering"), cex=1.4)
dev.off()
