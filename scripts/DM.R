##### Code for DM
##### 
##### Author: Jong Kyoung Kim (jkkim@ebi.ac.uk)
##### Last Update: 12/03/2015

##### Contents 
##### 1 Internal R functions
#####
##### 2 Examples

##### 1 Internal R functions

# or USE source("http://www.r-statistics.com/wp-content/uploads/2010/04/Quantile.loess_.r.txt")
Quantile.loess<- function(Y, X = NULL, 
                          number.of.splits = NULL,
                          window.size = 20,
                          percent.of.overlap.between.two.windows = NULL,
                          the.distance.between.each.window = NULL,
                          the.quant = .95,
                          window.alignment = c("center"), 
                          window.function = function(x) {quantile(x, the.quant)}, ...)
{
  if(!is.null(number.of.splits)) {window.size <- ceiling(length(Y)/number.of.splits)}
  if(is.null(the.distance.between.each.window)) {the.distance.between.each.window <- window.size}
  if(!is.null(percent.of.overlap.between.two.windows)) 
  {
    the.distance.between.each.window <- window.size * (1-percent.of.overlap.between.two.windows)
  }
  if(!require(zoo))   
  {
    print("zoo is not installed - please install it.")
    install.packages("zoo")
  }
  if(is.null(X)) {X <- index(Y)} 
  zoo.Y <- zoo(x = Y, order.by = X)
  new.Y <- rollapply(zoo.Y, width = window.size, 
                     FUN = window.function,
                     by = the.distance.between.each.window,
                     align = window.alignment)
  new.X <- attributes(new.Y)$index  
  new.Y.loess <- loess(new.Y~new.X, family = "sym",...)$fitted 
  return(list(y = new.Y, x = new.X, y.loess = new.Y.loess))
}

# DM for tag-based scRNA-seq data (e.g. UMI)
DM <- function(meanGenes, CV2Genes, windowSize=50) {
  meanGenesExpressed = meanGenes[meanGenes > 0 & !is.na(CV2Genes) & CV2Genes > 0]
  CV2GenesExpressed = CV2Genes[meanGenes > 0 & !is.na(CV2Genes) & CV2Genes > 0]
  qloess <- Quantile.loess(log10(CV2GenesExpressed), log10(meanGenesExpressed),
                           the.quant=.5, window.size=windowSize, percent.of.overlap.between.two.windows=.5)
  DM <- aaply(seq(1, length(meanGenesExpressed), 1), 1, function(x) {
    mindex<-which.min(abs(qloess[[2]]-log10(meanGenesExpressed[x])))
    as.vector(log10(CV2GenesExpressed[x]) - qloess[[1]][mindex])}, .expand=FALSE, .progress="text"
  )
  names(DM) <- names(meanGenesExpressed)
  DM
}


# 2 Examples
# 
# tag-based scRNA-seq data such as UMI
# Assume that you have a UMI count table (UMICountGenesSC2i), where row: Genes, column: Cells
# 
# meanGenes = rowMeans(UMICountGenesSC2i[rowMeans(UMICountGenesSC2i)>0,])
# CV2Genes = apply(UMICountGenesSC2i[rowMeans(UMICountGenesSC2i)>0,], 1, var) / meanGenes^2
# 
# DMCV2 = DM(meanGenes, CV2Genes) # DM for UMI data
# 
# read-based scRNA-seq data (Raw read counts (normalised by size factors), RPM, RPKM, or FPKM)
# Assume that you have a read count table (countGenes), where row: Genes, column: Cells
# library(DESeq)
# sizeFactor <- estimateSizeFactorsForMatrix(countGenes) # If you normalise the raw counts by size factors
# nCountGenes = t(t(countGenes) / sizeFactor)
# colnames(nCountGenes) = rep("OS25", ncol(nCountGenes))
# geneMean = rowMeans(nCountGenes)
# geneLength = read.table("GRCm38.74.gene_length.txt", sep="\t", stringsAsFactor=FALSE) # if you use transripts, it is the length of transcripts. If you use genes, it is the length of union of exons
# meanCutoff = 10
# DMOS25 = DMLengthAdjusted(nCountGenes, "OS25", geneLength, minCount=meanCutoff)
# 
