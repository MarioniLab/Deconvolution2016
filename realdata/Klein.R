# Normalizing Drop-Seq data with Summation strategy
#
# ---- Packages ----

library('DESeq2')
library('edgeR')
library('scran')
source("./functions.R") 

# ---- Data ----
esd0 <- read.csv("../data/GSM1599494_ES_d0_main.csv.bz2", header=FALSE, stringsAsFactors=FALSE, row.names=1)
esd7 <- read.csv("../data/GSM1599499_ES_d7_LIFminus.csv.bz2", header=FALSE, stringsAsFactors=FALSE, row.names=1)

colnames(esd0) <- paste0("d0.", seq_len(ncol(esd0)))
colnames(esd7) <- paste0("d7.", seq_len(ncol(esd7)))
counts <- merge (esd0, esd7, by=0, all=TRUE)
rownames(counts) <- counts[,1]
counts <- as.matrix(counts[,-1])

cell.num <- c(ncol(esd0), ncol(esd7))
timepoint <- rep(c("d0", "d7"), cell.num)
# timepointCol <- rep(c("red", "blue"), cell.num)

# ---- Gene-Filter ----
# Remove lowly expressed Genes

highE<- rowMeans(counts) >= 0.2
countsHE <- counts[highE,]
minExpr <- 1
featuresHE <- featureCalc(countsHE, htseq = FALSE, minExpr)

# ---- Size-Factors ----

#Normal SF
geoMeans_HE <- noZeroGM(countsHE)
szf_HE <- estimateSizeFactorsForMatrix(countsHE,geoMeans = geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHE),method="TMM") * colSums(countsHE)

#Normal Library Size
libfactor <- colSums(countsHE)

#Deconvoluted
szfCluster <- quickCluster(countsHE,deepSplit = 1)
szf_alClust <- computeSumFactors(countsHE,cluster=szfCluster)

# ---- SF-Difference ----

scaled_factors<- data.frame(DESeq=szf_HE/median(szf_HE), TMM=tmm/median(tmm), LibSize=libfactor/median(libfactor), Deconvolution=szf_alClust/median(szf_alClust), check.names=FALSE)
cell.col <- rgb(0,0,0) 
line.col <- "red"

cairo_pdf("Klein_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors,log="y",ylab="Size factors",cex.axis=1.5,cex.lab=1.8,xaxt='n')
text(labels=colnames(scaled_factors),xpd=TRUE,cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.1,srt=45, pos=2)
dev.off()

cairo_pdf("Klein_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="DESeq",ylab="Deconvolution",xlim=c(0.1,3.5),ylim=c(0.1,3.5),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

cairo_pdf("Klein_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$LibSize,scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="Library size",ylab="Deconvolution",xlim=c(0.1,8),ylim=c(0.1,8),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

cairo_pdf("Klein_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM,scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="TMM",ylab="Deconvolution",xlim=c(0.1,8),ylim=c(0.1,8),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

# ---- Saving for DE and other analyses ---

saveRDS(list(Cells=data.frame(Cell=colnames(counts), Time=timepoint), Counts=countsHE, SF=scaled_factors), file="KleinData.rds")

