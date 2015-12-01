CV2 <- function(x) {
    #Function that calculates the squared CV for a given x
  (sd(x)/mean(x))^2
}

featureCalc <- function(counts,htseq,exprmin,hvg) {
    #Takes countdata, the logical variables htseq (if True the htseq statistics are removed) and expression minium (usually 1) and a logical vector indicating whether a gene is highly variable
    #The output is a matrix with features as rownames and different properties of the features as coloumns
  if (htseq) {
    counts <- counts[ - grep("__",rownames(counts)),] ##remove htseq statistics
  }
  spike <- rep(FALSE,nrow(counts))
  spike[grep("ERCC",rownames(counts))] <- TRUE
  cv2 <- apply(counts,1,CV2)
  mean <- apply(counts,1,mean)
  expInCells <- apply(counts,1, function(x) length(x[x >= exprmin]))
  if (missing(hvg)) {
    hvg <- rep(NA,nrow(counts))
  }
  feats <- data.frame(spike,cv2,mean,hvg,expInCells,stringsAsFactors = F)
  rownames(feats) <- rownames(counts)
  return(feats)
}

drawBoxplot <- function(cnts,grps,gene) {
    #Takes in count matrix cnts, a vector containg the cluster number for each cell in SAME ORDER as the coutn matrix  and a gene name to plot as character
    # returns ggplot
    stopifnot(length(colnames(cnts)) == length (grps))
    stopifnot(length(which(rownames(cnts) == gene)) == 1)
    cnts <- data.frame(t(cnts),"cluster"=factor(grps))
    p <- ggplot(cnts, aes_string("cluster",gene))
    p <- p + geom_boxplot(aes(fill=cluster)) + geom_jitter(aes(colour=cluster), alpha = 0.5,size = 0.9) + theme(panel.background = element_blank(),legend.position="none")
}

ClusterSums <- function(cnts, clusters, sz) {
    ## Function that takes a count matrix, a vector containing cluster identites for coloumns of the cnts and a desired number (sz) of cells to sum up
    oldW <- getOption("warn")
    options(warn = -1)
    stopifnot(ncol(cnts) == length(clusters) & class(clusters) == "factor")
    nclust <- length(levels(clusters))
    Names <- levels(clusters)
    newcounts <- NULL
    identities <- NULL
    results <- list(newcounts,identities)
    for (i in c(1:nclust))  {
        Subset <- cnts[, clusters == Names[i]]
        ncells <- ncol(Subset)
        nsums <- ncells %/% sz
        nrest <- ncells %% sz
        shuffcells <- sample(colnames(Subset),ncol(Subset))
        cells2sum <- if(nrest == 0) split(shuffcells,rep(1:(ncol(Subset)/sz), each = sz)) else split(shuffcells,rep(1:((ncol(Subset)+(sz-nrest))/sz), each = sz))
        for (j in c(1:(nsums + (nrest != 0)))) {
            idnt <- cells2sum[[j]]
            sumcounts <- rowSums(Subset[,cells2sum[[j]],drop =FALSE])/length(idnt)
            results[[1]] <- cbind(results[[1]],sumcounts)
            colnames(results[[1]])[ncol(results[[1]])] <- paste("C",as.character(i),"S",as.character(j),sep="_") 
            results[[2]] <- c(results[[2]],list(idnt))
            names(results[[2]])[length(results[[2]])] <- paste("C",as.character(i),"S",as.character(j),sep="_") 
        }

    }
    rownames(results[[1]]) <- rownames(cnts)
    return(results)
    options(warn = oldW)
}

clusterForSummation<- function(counts, minSize = 200){
    ## This function generates a cluster vector containing the cluster number assigned to each cell. It takes the counts matrix and a minimum number of Cells per cluster as input. The minimum number should be at least twice as large as the largest group used for summation.
    require(dynamicTreeCut)
    stopifnot(ncol(counts) > minSize)
    if (ncol(counts) < 2 * minSize) {
        minSize <- as.integer((ncol(counts) / 5L))
        print(paste("Cluster size is scaled down to", minSize, ". You might need to adjust the 'sizes' argument for normalizeBySums accordingly."))
    }
    distM <- as.dist( 1 - cor(counts, method = 'spearman'))
    htree <- hclust(distM, method = 'ward.D2')
    clusters <- factor(unname(cutreeDynamic(htree, minClusterSize = minSize, method = 'hybrid', distM = as.matrix(distM), deepSplit = 0,pamStage = TRUE, verbose = 0, respectSmallClusters = TRUE)))
    if ( levels(clusters)[1] == 0) {
        clusters[clusters == 0] <- NA
        print(paste(sum(is.na(clusters)), "Cells could not be assigned to any cluster and are left unassaigned (NA). You might want to consider removing these cells from further analysis."))
    }
    return(clusters)
}
