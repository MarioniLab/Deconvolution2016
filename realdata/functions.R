CV2 <- function(x) {
    #Function that calculates the squared CV for a given x
  (sd(x)/mean(x))^2
}

noZeroGM <- function(x) {
    # Calculates the geometric mean after removing zeros from each row.
    x <- log(x)
    x[!is.finite(x)] <- NA_real_
    exp(rowMeans(x, na.rm=TRUE))
}

featureCalc <- function(counts,htseq,exprmin,hvg) {
    #Takes countdata, the logical variables htseq (if True the htseq statistics are removed) and expression minium (usually 1) and a logical vector indicating whether a gene is highly variable
    #The output is a matrix with features as rownames and different properties of the features as coloumns
  if (htseq) {
    counts <- counts[ - grep("__",rownames(counts)),] ##remove htseq statistics
  }
  spike <- grepl("ERCC",rownames(counts))
  cv2 <- apply(counts,1,CV2)
  mean <- rowMeans(counts)
  expInCells <- rowSums(counts >= exprmin)
  if (missing(hvg)) {
    hvg <- rep(NA,nrow(counts))
  }
  feats <- data.frame(spike,cv2,mean,hvg,expInCells,stringsAsFactors = FALSE)
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

compareHVG <- function(HVGdf,topX) {
    #Takes data.frame with colomns represented orderd gene names by highes DM value and an integer for the top X genes..
    #Returns matrix showing overlap between each method
    HVGdf <- HVGdf[c(1:topX),]
    result <- NULL 
    for (i in c(1:ncol(HVGdf))) {
             overlap <- sapply(HVGdf,function(x) length(intersect(HVGdf[,i],x))) ## Calc number of overlap between HVG[i] and all other HVGs
             uniq <- length(HVGdf[,i]) - length(intersect(HVGdf[,i],unlist(HVGdf[,-i]))) ## Calculate number of unique HVGs
             result <- cbind(result,c(overlap,uniq))
    }
    colnames(result) <- colnames(HVGdf)
    rownames(result) <- c(colnames(HVGdf),"Unique")
    return(result)
}
