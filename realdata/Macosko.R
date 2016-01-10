# scran Normalization Macosko

# ---- Packages ----
pkgs <- "/nfs/research2/marioni/Karsten/utils/Rpack/"
library('rJava',lib.loc=pkgs)
library('xlsxjars',lib.loc=pkgs)
library('xlsx',lib.loc=pkgs)
library('VennDiagram',lib.loc=pkgs)
library('zoo',lib.loc=pkgs)
library('tsne',lib.loc=pkgs)
library('DESeq2')
library('plyr')
library('edgeR')
library('dynamicTreeCut',lib.loc=pkgs)
library('moduleColor',lib.loc=pkgs)
library('lattice',lib.loc=pkgs)
library('ggplot2',lib.loc=pkgs)
library('gridExtra',lib.loc=pkgs)
library('pheatmap',lib.loc=pkgs)
library("RColorBrewer")
library('scran')
source("./functions.R") 
source("./DM.R")

#load data (takes some time)
Retina1 <- read.delim("../data/GSM1626793_P14Retina_1.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina1) <- c(1:ncol(Retina1))
Retina2 <- read.delim("../data/GSM1626794_P14Retina_2.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina2) <- c(1:ncol(Retina2)) + ncol(Retina1)
Retina3 <- read.delim("../data/GSM1626795_P14Retina_3.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina3) <- c(1:ncol(Retina3)) + ncol(Retina1) + ncol(Retina2)
Retina4 <- read.delim("../data/GSM1626796_P14Retina_4.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina4) <- c(1:ncol(Retina4)) + ncol(Retina1) + ncol(Retina2) + ncol(Retina3)
Retina5 <- read.delim("../data/GSM1626797_P14Retina_5.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina5) <- c(1:ncol(Retina5)) + ncol(Retina1) + ncol(Retina2) + ncol(Retina3) + ncol(Retina4)
Retina6 <- read.delim("../data/GSM1626798_P14Retina_6.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina6) <- c(1:ncol(Retina6)) + ncol(Retina1) + ncol(Retina2) + ncol(Retina3) + ncol(Retina4) + ncol(Retina5)
Retina7 <- read.delim("../data/GSM1626799_P14Retina_7.digital_expression.txt", stringsAsFactors=FALSE)
colnames(Retina7) <- c(1:ncol(Retina7)) + ncol(Retina1) + ncol(Retina2) + ncol(Retina3) + ncol(Retina4) + ncol(Retina5) + ncol(Retina6)

retina <- merge(Retina1, Retina2, by = 1, all=T)
retina <- merge(retina, Retina3, by = 1, all=T)
retina <- merge(retina, Retina4, by = 1, all=T)
retina <- merge(retina, Retina5, by = 1, all=T)
retina <- merge(retina, Retina6, by = 1, all=T)
retina <- merge(retina, Retina7, by = 1, all=T)
rownames(retina) <- retina[,1]
retina <- retina[,2:ncol(retina)]
retina <- na.omit(retina)
rm("Retina1","Retina2","Retina3","Retina4","Retina5","Retina6","Retina7")

# Only choose training set as definied by authors which removes cells with very low counts
exprGenes <- apply(retina,2,function(x) sum(x != 0))
subst <- exprGenes >= 900
retinaSs <- retina[,subst]


set.seed(1)
smpl <- sample(1:ncol(retinaSs),4000,replace = FALSE)
retina.smpl <- retinaSs[,smpl]

# ---- Gene-Filter ----

retina.smpl <- na.omit(retina.smpl)

IncCrit <- ncol(retina.smpl)/20

highE<- rowSums(retina.smpl) >= IncCrit
countsHE <- retina.smpl[highE,]

# ---- Size-Factors ----

#Normal SF
geoMeans_HE <- apply(countsHE, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
szf_HE <- estimateSizeFactorsForMatrix(countsHE,geoMeans = geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHE),method="TMM") * colSums(countsHE)

#Normal Library Size
libfactor <- colSums(countsHE) / 10000

#Deconvoluted
szfCluster <- quickCluster(countsHE)
szf_alClust <- normalizeBySums(countsHE,cluster=szfCluster,positive = TRUE)

# ---- SF-Difference ----
cellCols <- labels2colors(as.integer(szfCluster))
scaled_factors<- data.frame("DESeq" = szf_HE / median(szf_HE),"TMM" = tmm/median(tmm), "Library size" = libfactor/median(libfactor), "Deconvolution" = szf_alClust/median(szf_alClust),check.names=FALSE)

cairo_pdf("Macosko_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors,log="y",ylab="Size factors",cex.axis=1.5,cex.lab=1.8,xaxt='n')
text(labels=colnames(scaled_factors),xpd=TRUE,cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.1,srt=45, pos=2)
dev.off()

cairo_pdf("Macosko_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,scaled_factors$Deconvolution,cex=0.6,pch=19,col="#00000073",xlab="DESeq",ylab="Deconvolution",xlim=c(0.1,3.5),ylim=c(0.1,20),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")
dev.off()

cairo_pdf("Macosko_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors[,"Library size"],scaled_factors$Deconvolution,cex=0.6,pch=19,col=cellCols,xlab="Library size",ylab="Deconvolution",xlim=c(0.1,20),ylim=c(0.1,20),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")
dev.off()

cairo_pdf("Macosko_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors[,"TMM"],scaled_factors$Deconvolution,cex=0.6,pch=19,col="#00000073",xlab="TMM",ylab="Deconvolution",xlim=c(0.1,20),ylim=c(0.1,20),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")
dev.off()

distM <- as.dist(1-cor(countsHE,method="spearman"))
tsn <- tsne(distM)
plot(tsn,col=cellCols)

counts.pca <- prcomp(t(log(countsHE+1)))
plot(counts.pca$x[,1:2],col=cellCols,pch=19,cex=0.6)
tsn <- tsne(counts.pca$x[,1:10],initial_dims = 10,k=2)
cairo_pdf("TSNE-Plot.pdf")
plot(tsn,col=cellCols,pch=19,cex=0.6)
dev.off()
