# scran Normalization

# ---- Packages ----
library('DESeq2')
library('edgeR')
library('scran')
source("./functions.R") 
#source("./DM.R")

# ---- Data ----

readFormat <- function(infile) { 
    metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] # First column is empty.
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
    metadata <- as.data.frame(t(metadata))
    counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] # First column after row names is some useless filler.
    counts <- as.matrix(counts)
    return(list(metadata=metadata, counts=counts))
}

cellular <- readFormat("../data/expression_mRNA_17-Aug-2014.txt")
erccspikes <- readFormat("../data/expression_spikes_17-Aug-2014.txt")
mitochondrial <- readFormat("../data/expression_mito_17-Aug-2014.txt")

stopifnot(identical(cellular$metadata$cell_id, erccspikes$metadata$cell_id)) # should be the same.
m <- match(cellular$metadata$cell_id, mitochondrial$metadata$cell_id)
mitochondrial$metadata <- mitochondrial$metadata[m,]
mitochondrial$counts <- mitochondrial$counts[,m]
stopifnot(all(cellular$metadata$cell_id==mitochondrial$metadata$cell_id)) # should now be the same.

# Combine data sets into one expression matrix 
counts <- rbind(cellular$counts,mitochondrial$counts,erccspikes$counts)
metadata <- cellular$metadata
minExpr <- 1

features <- featureCalc(counts = counts,htseq = FALSE, exprmin = minExpr) # calcuate feature data 
features$mito <- grepl("mt-",rownames(counts)) # Add logical coloumn with mt gene indices

# ---- Gene-Filter ----

# Remove lowly expressed Genes

highE<- rowMeans(counts) >= 0.2
countsHE <- counts[highE,]

featuresHE <- featureCalc(countsHE, htseq = FALSE, minExpr)

# ---- Size-Factors ----

countsHEnoSpike <- countsHE[!featuresHE$spike,]

#Normal SF
geoMeans_HE <- noZeroGM(countsHEnoSpike)
szf_HE <- estimateSizeFactorsForMatrix(countsHEnoSpike, geoMeans = geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHEnoSpike),method="TMM") * colSums(countsHEnoSpike)

#Normal Library Size
libfactor <- colSums(countsHEnoSpike) 

#ERCC Size factors
erccfactor <- colSums(countsHE[featuresHE$spike,])

#Deconvoluted
szfCluster <- quickCluster(countsHEnoSpike)
szf_alClust <- normalizeBySums(countsHEnoSpike, clusters = szfCluster)

# ---- SF-Difference ----

scaled_factors <- data.frame(DESeq=szf_HE/median(szf_HE), TMM=tmm/median(tmm), LibSize=libfactor/median(libfactor), Deconvolution=szf_alClust/median(szf_alClust), 
                             Spikes=erccfactor/median(erccfactor), check.names=FALSE)

#Comparison of SF Distribution
#Colors
line.col <- "red"
olig.col <- rgb(1,0.5,0)
pyr.col <- "dodgerblue"
cell.col <- rgb(0,0,0) 

setCols <- rep(cell.col, nrow(metadata))
setCols[metadata$level1class=="oligodendrocytes"] <- olig.col
setCols[metadata$level1class=="pyramidal CA1"] <- pyr.col

set.seed(1)
shffl <- sample(ncol(countsHE)) # Shuffling cells to avoid one colour plotting over another.

cairo_pdf("Zeisel_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors[,-ncol(scaled_factors)],log="y",ylab="Size factors",cex.axis=1.5,cex.lab=1.8,xaxt='n')
text(labels=colnames(scaled_factors[,-ncol(scaled_factors)]),xpd=TRUE,cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.1,srt=45, pos=2)
dev.off()

#SF vs Deconvoluted

cairo_pdf("Zeisel_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq[shffl],scaled_factors$Deconvolution[shffl],xlim=c(0.1,6),ylim=c(0.1,6),pch=16,xlab="DESeq",ylab="Deconvolution",log="xy",cex.axis=1.5,cex.lab=1.8,col=setCols[shffl])
abline(0,1,col=line.col,lwd=2)
dev.off()

#Libfactor vs Deconvoluted
cairo_pdf("Zeisel_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$LibSize[shffl],scaled_factors$Deconvolution[shffl],xlim=c(0.1,7),ylim=c(0.1,7),pch=16,xlab="Library size",ylab="Deconvolution",log="xy",cex.axis=1.5,cex.lab=1.8,col=setCols[shffl])
abline(0,1,col=line.col,lwd=2)
dev.off()

#TMM vs Deconvoluted
cairo_pdf("Zeisel_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM[shffl],scaled_factors$Deconvolution[shffl],xlim=c(0.1,7),ylim=c(0.1,7),pch=16,col=setCols[shffl],xlab="TMM",ylab="Deconvolution",cex.axis=1.5,cex.lab=1.8,log="xy")
abline(0,1,col=line.col,lwd=2)
dev.off()

# ---- Spike-Ins ----

cairo_pdf("Zeisel_ERCCvSF.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,scaled_factors$Spikes,ylab="Spike-in",xlab="DESeq",log="xy",pch=16,col="#000000",cex.axis=1.5,cex.lab=1.8)
dev.off()
cairo_pdf("Zeisel_ERCCvTMM.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM,scaled_factors$Spikes,ylab="Spike-in",xlab="TMM",log="xy",cex=0.6,pch=19,col="#000000",cex.axis=1.5,cex.lab=1.8)
dev.off()
cairo_pdf("Zeisel_ERCCvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$Deconvolution,scaled_factors$Spikes,ylab="Spike-in",xlab="Deconvolution",log="xy",cex=0.6,pch=19,col="#000000",cex.axis=1.5,cex.lab=1.8)
dev.off()
cairo_pdf("Zeisel_ERCCvLib.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$LibSize,scaled_factors$Spikes,ylab="Spike-in",xlab="Library size",log="xy",cex=0.6,pch=19,col="#000000",cex.axis=1.5,cex.lab=1.8)
dev.off()

# ---- Saving for DE and other analyses ---

saveRDS(list(Cells=metadata, Counts=countsHEnoSpike, SF=scaled_factors), file="ZeiselData.rds")
