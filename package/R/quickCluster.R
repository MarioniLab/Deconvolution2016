quickCluster <- function(counts, minSize = 200, ...)
# This function generates a cluster vector containing the cluster number assigned to each cell.
# It takes the counts matrix and a minimum number of Cells per cluster as input.
# The minimum number should be at least twice as large as the largest group used for summation.
#
# written by Karsten Bach
# created 1 December 2015
{   
    if (ncol(counts) <=  minSize){
        stop('Less Cells than mininimum cluster size')
    }
    if (ncol(counts) < 2 * minSize) {
        minSize <- as.integer((ncol(counts) / 5L))
        warning(paste("MinSize scaled down to", minSize))
    }
    distM <- as.dist( 1 - cor(counts, method = 'spearman'))
    htree <- hclust(distM, method = 'ward.D2')
    clusters <- factor(unname(cutreeDynamic(htree, minClusterSize = minSize,distM = as.matrix(distM),verbose = 0, ...)))
    if ( levels(clusters)[1] == 0) {
        clusters[clusters == 0] <- NA
        warning(paste(sum(is.na(clusters)), "cells are left unassaigned (NA)"))
    }
    return(clusters)
}
