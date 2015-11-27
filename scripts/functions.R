CV2 <- function(x) {
    #Function that calculates the squared CV for a given x
  (sd(x)/mean(x))^2
}

percentofreads <- function(dta,ssdta) {
    ## Function used in cellCalc to determine the percentage of reads that are within a subset of the features (e.g only Exons, ERCC etc.) 
  if (dim(dta)[2] != dim(ssdta)[2]  ){
    print("error")
    break
  }
  sums_total <- unname(colSums(dta))
  sums_ss_data <- unname(colSums(ssdta))
  read_fraction <- sums_ss_data / sums_total
  return(read_fraction)
}

cellCalc <- function (counts,mitogenes,readcounts,exprmin,celltype) {
  stopifnot(class(counts) == "data.frame", nrow(readcounts) == ncol(counts), length(celltype) == ncol(counts), length(which(colnames(counts) %in% readcounts[,1])) == ncol(counts))
  ##Creates CellData matrix containing various QC information on each cell
  ## Takes countData, List of mitochondrial Genes, # of mapped reads and celltype vector for each cell and an exprmin
  ERCC_ss_data <- counts[grep("ERCC-",rownames(counts)),]
  Exons_ss_data <- counts[grep("ENSMUSG",rownames(counts)),]
  not_aligned_ss_data <- counts[grep("__not_aligned",rownames(counts)),]
  lowQ_ss_data <- counts[grep("__too_low_aQual",rownames(counts)),]
  ambiguous_ss_data <- counts[grep("__ambiguous",rownames(counts)),]
  no_feature_ss_data <- counts[grep("__no_feature",rownames(counts)),]
  mito_ss_data <-counts[rownames(counts) %in% mitogenes,]
  ERCC_prcnt <- percentofreads(counts,ERCC_ss_data)
  Exons_prcnt <- percentofreads(counts,Exons_ss_data)
  not_aligned_prcnt <- percentofreads(counts,not_aligned_ss_data)
  lowQ_prcnt <- percentofreads(counts,lowQ_ss_data)
  ambiguous_prcnt <-percentofreads(counts,ambiguous_ss_data)
  no_feature_prcnt <- percentofreads(counts,no_feature_ss_data)
  total_UMI <- colSums(counts)
  total_UMI_genes <- colSums(Exons_ss_data)
  mito_prcnt_of_genes <- percentofreads(counts,mito_ss_data)
  selection_vector <- c(1:dim(counts)[2])
  expr_genes <- apply(Exons_ss_data,2, function(x) length(x[x >= exprmin]))
  QC_set <- data.frame(ERCC_prcnt,not_aligned_prcnt,lowQ_prcnt,ambiguous_prcnt,no_feature_prcnt,Exons_prcnt,total_UMI,mito_prcnt_of_genes,total_UMI_genes,selection_vector,expr_genes,celltype,stringsAsFactors = F)
  rownames(QC_set) <- colnames(counts)
  QC_set <- merge(QC_set,readcounts,by.x = "row.names" , by.y = "names")
  rownames(QC_set) <- QC_set$Row.names
  QC_set <- QC_set[,-1]
  return(QC_set)
}

cellCalc_nonhtseq <- function(counts,celltype,exprmin) {
  ##Creates CellData matrix containing various QC information on each cell
  ## Takes countData, List of mitochondrial Genes, # of mapped reads and celltype vector for each cell and an exprmin
  stopifnot(class(counts) == "data.frame", length(celltype) == ncol(counts))
  total_UMI <- colSums(counts)
  expr_genes <- apply(counts,2,function(x) length(x[x >= exprmin]))
  QC_set <- data.frame(celltype,total_UMI,expr_genes,stringsAsFactors = F)
  rownames(QC_set) <- colnames(counts)
  return(QC_set)
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
    


  QCplot <- function (QCset) {
     ## function that takes a QCset produced with cellCalc as input and gives a QC plot as output
     ## Only applicable to data that has been produced by htseq-count 
    oldW <- getOption("warn")
    options(warn = -1)
    old_par <- par()
    QCsetOrd <- QCset[order(QCset$Exons_prcnt),] 
    par(mar=c(5.1,4.1,4.1,4))
    barplot(QCsetOrd$ERCC_prcnt+QCsetOrd$not_aligned_prcnt+QCsetOrd$lowQ_prcnt+QCsetOrd$ambiguous_prcnt+QCsetOrd$no_feature_prcnt+QCsetOrd$Exons_prcnt,col='black',ylab='Percent of all UMIs',xlab='Cells',ylim=c(0,1),xlim=c(-1,nrow(QCset)+12),width=1,space=0)
    par(new=TRUE)
    barplot(QCsetOrd$ERCC_prcnt+QCsetOrd$not_aligned_prcnt+QCsetOrd$ambiguous_prcnt+QCsetOrd$no_feature_prcnt+QCsetOrd$Exons_prcnt,col='powderblue',axes=FALSE,ylim=c(0,1),xlim=c(-1,nrow(QCset)+12),width=1,space=0)
    par(new=TRUE)
    barplot(QCsetOrd$ERCC_prcnt+QCsetOrd$not_aligned_prcnt+QCsetOrd$no_feature_prcnt+QCsetOrd$Exons_prcnt,col='darkorchid',axes=FALSE,ylim=c(0,1),xlim=c(-1,nrow(QCset)+12),width=1,space=0)
    par(new=TRUE)
    barplot(QCsetOrd$ERCC_prcnt+QCsetOrd$no_feature_prcnt+QCsetOrd$Exons_prcnt,col='salmon',axes=FALSE,ylim=c(0,1),xlim=c(-1,nrow(QCset)+12),width=1,space=0)
    par(new=TRUE)
    barplot(QCsetOrd$ERCC_prcnt+QCsetOrd$Exons_prcnt,col='deepskyblue',axes=FALSE,ylim=c(0,1),xlim=c(-1,nrow(QCset)+12),width=1,space=0)
    par(new=TRUE)
    barplot(QCsetOrd$Exons_prcnt,col='seagreen1',axes=FALSE,ylim=c(0,1),xlim=c(-1,nrow(QCset)+12),width=1,space=0)
    par(new=TRUE)
    plot(1:nrow(QCsetOrd),log10(QCsetOrd$readcounts),axes = FALSE, pch = 16, xlab ="", ylab = "", xlim = c(par("usr")[1],par("usr")[2]), ylim=c(log10(min(QCsetOrd$readcounts)),log10(max(QCsetOrd$readcounts))),xaxs="i")
    axis(side=4)
    axis(side=1)
    legend("bottom",inset =c(-0.15,-0.2),x.intersp=0.5,y.intersp=0.5,horiz=T,xpd = T,bty="n",cex=0.9,pch=c(15,15,15,15,15,15,16),col=c("seagreen1","deepskyblue","salmon","darkorchid","powderblue","black","black"), legend = c("Exons","ERCC","no_feat","not_aligned","ambiguous","lowQ","reads"))
    mtext("Log10 readcounts",side=4,line=3)
    par(old_par)
    options(warn = oldW)
  }

plotgroups <- function(hier1,Color,cellmeta) 
{
  if (length(hier1$order) != length(Color) | length(hier1$order) != length(cellmeta) ) 
        { 
         stop("ERROR: length of color vector not compatible with no. of objects in the hierarchical tree");
    } else {
                    barplot(height=rep(-0.05*max(hier1$height), length(Color)), col= as.character(Color[hier1$order]),
                    border = NA, main="",space=0,offset = - 3, axes = FALSE,add = TRUE)
                    barplot(height=rep(-0.05*max(hier1$height), length(cellmeta)), col=labels2colors(as.numeric(cellmeta[hier1$order])),
                    border=NA, main= "", space = 0, offset = - 4, axes = FALSE, add = TRUE)
        
      }
}


compareClust <- function(grps1,grps2) {
    # Function that takes two vectors with group names for cells, the grps1[i] represents the group and i the cell number according to the original sorting in the matrix. !! Important that the two vectors have the same order, that is the one of the original data.frame object.
    # Output is a matrix with coloumn representing Clusters1 and rows representing Clusters2
    stopifnot(length(grps1) == length(grps2))
    grps1 <- factor(grps1)
    grps2 <- factor(grps2)
    lvls1 <- levels(grps1)
    lvls2 <- levels(grps2)
    result <- NULL 
    for (i in c(1:length(lvls1))) {
             overlap <- NULL
             members1 <- which(grps1 == lvls1[i])
             for (j in c(1:length(lvls2))) {

                members2 <- which(grps2 == lvls2[j])
                overlap[j] <- sum(members1 %in% members2) / length(union(members1,members2))
             }
             result <- cbind(result,overlap)
    }
    colnames(result) <- lvls1
    rownames(result) <- lvls2
    result
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

cellCor <- function(count1,count2) {
    # simple function to calculate colwise correlation coefficients
    stopifnot(ncol(count1) == ncol(count2) & nrow(count1) == nrow(count2))
    correlation <- 0
        for (i in 1:ncol(count1)) {
           correlation[i] <- cor(count1[,i],count2[,i])
        }
    correlation
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

