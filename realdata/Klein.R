# Normalizing Drop-Seq data with Summation strategy
#
# ---- Packages ----
pkgs <- "/nfs/research2/marioni/Karsten/utils/Rpack/"
library('rJava',lib.loc=pkgs)
library('xlsxjars',lib.loc=pkgs)
library('xlsx',lib.loc=pkgs)
library('zoo',lib.loc=pkgs)
library('DESeq2')
library('plyr')
library('edgeR')
library('dynamicTreeCut',lib.loc=pkgs)
library('lattice',lib.loc=pkgs)
library('ggplot2',lib.loc=pkgs)
library('gridExtra',lib.loc=pkgs)
library('pheatmap',lib.loc=pkgs)
library("RColorBrewer")
library('AK47')
source("./functions.R") 
source("./DM.R")

# ---- Data ----
esd0 <- read.csv("../data/GSM1599494_ES_d0_main.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  
esd7 <- read.csv("../data/GSM1599499_ES_d7_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  

colnames(esd0) <- c("feature",c(1:933))
colnames(esd7) <- c("feature",c(1920:2717))

counts <- merge (esd0,esd7,by="feature",all=T)
rownames(counts) <- counts$feature
counts <- counts[,!(colnames(counts) %in% "feature")]

celltypeCol <- c(rep("red",933),rep("blue",798))
celltype <- c(rep("d0",933),rep("d7",798))

# ---- Gene-Filter ----
# Remove lowly expressed Genes
IncCrit <- ncol(counts)/5

highE<- rowSums(counts) >= IncCrit
countsHE <- counts[highE,]
minExpr <- 1
featuresHE <- featureCalc(countsHE, htseq = FALSE, minExpr)


# ---- Size-Factors ----

#Normal SF
geoMeans_HE <- apply(countsHE, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
szf_HE <- estimateSizeFactorsForMatrix(countsHE,geoMeans = geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHE),method="TMM") * colSums(countsHE)

#Normal Library Size
libfactor <- colSums(countsHE) / 10000

#Deconvoluted
szfCluster <- quickCluster(countsHE,deepSplit = 1)
szf_alClust <- normalizeBySums(countsHE,cluster=szfCluster)


# ---- SF-Difference ----
scaled_factors<- data.frame("DESeq" = szf_HE / median(szf_HE),"TMM" = tmm/median(tmm), "Library size" = libfactor/median(libfactor), "Deconvolution" = szf_alClust/median(szf_alClust),check.names=FALSE)
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
plot(scaled_factors[,"Library size"],scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="Library size",ylab="Deconvolution",xlim=c(0.1,8),ylim=c(0.1,8),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

cairo_pdf("Klein_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors[,"TMM"],scaled_factors$Deconvolution,pch=16,col=cell.col,xlab="TMM",ylab="Deconvolution",xlim=c(0.1,8),ylim=c(0.1,8),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col=line.col)
dev.off()

# ---- Normalization ----

counts.sf <- as.data.frame(t(t(countsHE) / szf_HE))
counts.tmm <- as.data.frame(t(t(countsHE) / tmm))
counts.lib <- as.data.frame(t(t(countsHE) / libfactor))
counts.decon <- as.data.frame(t(t(countsHE) / szf_alClust))

features.sf <- featureCalc(counts.sf, htseq=FALSE,minExpr)
features.tmm <- featureCalc(counts.tmm, htseq=FALSE,minExpr)
features.lib <- featureCalc(counts.lib, htseq=FALSE,minExpr)
features.decon <- featureCalc(counts.decon, htseq=FALSE,minExpr)

# ---- HVG-calling ----

dm.sf <- DM(meanGenes=features.sf$mean,CV2Genes=features.sf$cv2)
dm.tmm <- DM(meanGenes=features.tmm$mean,CV2Genes=features.tmm$cv2)
dm.lib <- DM(meanGenes=features.lib$mean,CV2Genes=features.lib$cv2)
dm.decon <- DM(meanGenes=features.decon$mean,CV2Genes=features.decon$cv2)

features.sfOrd <- features.sf[order(-dm.sf),]
features.tmmOrd <- features.tmm[order(-dm.tmm),]
features.libOrd <- features.lib[order(-dm.lib),]
features.deconOrd <- features.decon[order(-dm.decon),]

for (x in c(500,1000,2000)) {

    HVGlist <- list("DESeq"=rownames(features.sfOrd)[1:x], "TMM"=rownames(features.tmmOrd)[1:x], "Lib"=rownames(features.libOrd)[1:x], "Deconvolution"=rownames(features.deconOrd)[1:x])

    comparisonMatrix <- compareHVG(HVGlist)
    write.table(comparisonMatrix,paste("HVGcomparisonKlein",x,".txt",sep=""))
}
