# Normalizing Kowalczyk data with Summation strategy
#
# ---- Packages ----

require('DESeq2')
require('edgeR')
require('AK47')
require('openxlsx')
source("./functions.R") 
source("./DM.R")

# ---- Data ----

#load Black6
counts.log <- read.xlsx("../data/GSE59114_C57BL6_GEO_all.xlsx", startRow=2)
counts.log <- counts.log[, -c(2,1431:1436)] #removes UCSC transcripts and populations
counts.log$Gene.Symbol <- gsub("'", "", counts.log$Gene.Symbol)

rownames(counts.log) <- make.names(counts.log$Gene.Symbol, unique=TRUE)
counts.log <- counts.log[,-1]
counts <-2^counts.log - 1

# ---- Gene-Filter ----

highE <- rowMeans(counts) >= 0.2
countsHE <- as.matrix(counts[highE,]) 


# ---- Size-Factors ----

#Normal SF
tmp <- log(countsHE)
tmp[!is.finite(tmp)] <- NA_real_
geoMeans_HE <- exp(rowMeans(tmp, na.rm=TRUE))
szf_HE <- estimateSizeFactorsForMatrix(countsHE,geoMeans = geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHE),method="TMM") * colSums(countsHE)

#Normal Library Size
libfactor <- colSums(countsHE)

#Deconvoluted
szfCluster <- quickCluster(countsHE,deepSplit = 1)
szf_alClust <- normalizeBySums(countsHE,cluster=szfCluster)

# ---- SF-Difference ----

scaled_factors<- data.frame(DESeq=szf_HE/median(szf_HE), TMM=tmm/median(tmm), LibSize=libfactor/median(libfactor), Deconvolution=szf_alClust/median(szf_alClust), check.names=FALSE)
cell.col <- rgb(0,0,0) 
line.col <- "red"

cairo_pdf("Kowalczyk_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors,log="y",ylab="Size factors",cex.axis=1.5,cex.lab=1.8,xaxt='n')
text(labels=colnames(scaled_factors),xpd=TRUE,cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.22,srt=45, pos=2)
dev.off()

cairo_pdf("Kowalczyk_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="DESeq",ylab="Deconvolution",xlim=c(0.1,2),ylim=c(0.1,2),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

cairo_pdf("Kowalczyk_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$LibSize,scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="Library size",ylab="Deconvolution",xlim=c(0.1,2),ylim=c(0.1,2),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

cairo_pdf("Kowalczyk_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM,scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="TMM",ylab="Deconvolution",xlim=c(0.1,2),ylim=c(0.1,2),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

# ---- Saving for DE and other analyses ---

saveRDS(list(Cell=colnames(counts), Counts=countsHE, SF=scaled_factors), file="KowalczykData.rds")


