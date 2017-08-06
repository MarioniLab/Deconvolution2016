generateRawMeans <- function(ngenes=10000, ncells=200, mode=c("UMI", "read")) {
    mode <- match.arg(mode)
    if (mode=="UMI") {
        # Guestimated from inDrop data.
        true.means <- rgamma(ngenes, 2, 2)
        dispersions <- rep(0.1, ngenes)
    } else {
        # Guestimated from ESpresso data (see PlateEffect supplementaries).
        true.means <- 2^runif(ngenes, 0, 12)
        dispersions <- 100/true.means^0.6
    }
    
    all.facs <- 2^rnorm(ncells, sd=1)
    return(list(fitted=outer(true.means, all.facs), 
                mean=true.means, dispersion=dispersions, sf=all.facs))
}

sampleCounts <- function(means, dispersions) {
    matrix(rnbinom(length(means), mu=means, size=1/dispersions), ncol=ncol(means))
}

makeSFPlot <- function(sf, truth, is.de=NULL, main="", col="black") {
    logfold <- log2(sf) - log2(truth)

    if (is.null(is.de)) {
        fitted <- MASS::rlm(logfold ~ 1, na.action=na.exclude)
        resids <- residuals(fitted)
    } else {
        fitted <- MASS::rlm(logfold[-is.de] ~ 1, na.action=na.exclude)
        resids <- logfold - coef(fitted)
    }

    all.range <- range(truth)
    col <- rep(col, length.out=length(sf))
    shuffle <- sample(length(sf))
    
    plot(truth[shuffle], (truth * 2^resids)[shuffle], xlim=all.range, ylim=all.range, 
         ylab="Estimated factors", xlab="True factors", log="xy", pch=16, 
         col=col[shuffle], cex.axis=1.5, cex.lab=1.8, main=main, cex.main=1.8)
    abline(0, 1, col="red")

    DE.err <- 2^sqrt(mean(resids[is.de]^2))-1
    err.formatted <- format(round(DE.err*100, 1), nsmall=1)
    legend("topleft", bty="n", cex=1.2,
           legend=paste0("DE error = ", err.formatted, "%"))
    return(DE.err)
}

runAllMethods <- function(counts) {
    require(edgeR)
    require(DESeq2)
    require(scran)

    # TMM with raw counts:
    tmm.sf <- calcNormFactors(counts) * colSums(counts)

    # TMM with averaged counts:
    combined <- cbind(rowSums(counts), counts)
    tmm2.sf <- calcNormFactors(combined) * colSums(combined)
    tmm2.sf <- tmm2.sf[-1]

    # Size factors with raw counts (must counter zeros for both geometric mean and for each library):
    logvals <- log(counts)
    logvals[is.infinite(logvals)] <- NA_real_
    gm <- exp(rowMeans(logvals, na.rm=TRUE))
    size.sf <- estimateSizeFactorsForMatrix(counts, geoMeans=gm)

    # Size factors with averaged counts (still removes zeros in each library):
    size2.sf <- estimateSizeFactorsForMatrix(counts, geoMeans=rowMeans(counts))

    # Size factors with an added pseudo-count.
    sizeP.sf <- estimateSizeFactorsForMatrix(counts+1)
                                
    # Size factors with a library size-adjusted pseudo-count.
    lib.size <- colSums(counts)
    pcounts <- t(t(counts) + lib.size/mean(lib.size))
    sizeP2.sf <- estimateSizeFactorsForMatrix(pcounts)

    # Library size normalization.
    lib.sf <- colSums(counts)

    # Size factors with summation, no abundance filtering for simplicity.
    final.sf <- computeSumFactors(counts, clusters=NULL, min.mean=0)

    # Size factors with clustering prior to summation:
    if (ncol(counts) >= 200) {
        emp.clusters <- quickCluster(counts)
    } else {
        emp.clusters <- rep(1, ncol(counts))
    }
    final2.sf <- computeSumFactors(counts, clusters=emp.clusters, min.mean=0)

#    # knn-based method
#    final3.sf <- computeSumFactors(counts, mode="experimental", min.mean=0)

    # Reporting all methods.
    return(list(TMM=tmm.sf, TMM.ave=tmm2.sf,
                DESeq.geo=size.sf, DESeq.ave=size2.sf, DESeq.pseudo=sizeP.sf, DESeq.pseudo.lib=sizeP2.sf,
                Lib=lib.sf,
                Deconv=final.sf, Deconv.clust=final2.sf)) #Deconv.kNN=final3.sf))
}

