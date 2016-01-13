# scran Normalization Macosko

# ---- Packages ----

library('DESeq2')
library('edgeR')
library('scran')
source("./functions.R") 

#load data (takes some time)

Retina1 <- read.delim("../data/GSM1626793_P14Retina_1.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)
Retina2 <- read.delim("../data/GSM1626794_P14Retina_2.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)
Retina3 <- read.delim("../data/GSM1626795_P14Retina_3.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)
Retina4 <- read.delim("../data/GSM1626796_P14Retina_4.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)
Retina5 <- read.delim("../data/GSM1626797_P14Retina_5.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)
Retina6 <- read.delim("../data/GSM1626798_P14Retina_6.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)
Retina7 <- read.delim("../data/GSM1626799_P14Retina_7.digital_expression.txt.gz", skip=1, row.names=1, header=FALSE)

present.in.all <- Reduce(intersect, list(rownames(Retina1), rownames(Retina2), rownames(Retina3), rownames(Retina4), 
		                         rownames(Retina5), rownames(Retina6), rownames(Retina7)))
retina <- cbind(Retina1[present.in.all,], Retina2[present.in.all,], Retina3[present.in.all,], Retina4[present.in.all,], 
                Retina5[present.in.all,], Retina6[present.in.all,], Retina7[present.in.all,])
retina <- as.matrix(retina)
rm("Retina1","Retina2","Retina3","Retina4","Retina5","Retina6","Retina7")
gc()

# Only choose training set as definied by authors which removes cells with very low counts

exprGenes <- colSums(retina!=0)
subst <- exprGenes >= 900
retinaSs <- retina[,subst]

set.seed(1)
smpl <- sample(ncol(retinaSs), 4000, replace = FALSE)
retina.smpl <- retinaSs[,smpl]

# ---- Gene-Filter ----

highE<- rowSums(retina.smpl) >= ncol(retina.smpl)/20
countsHE <- retina.smpl[highE,]

# ---- Size-Factors ----

#Normal SF
geoMeans_HE <- noZeroGM(countsHE)
szf_HE <- estimateSizeFactorsForMatrix(countsHE, geoMeans=geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHE), method="TMM") * colSums(countsHE)

#Normal Library Size
libfactor <- colSums(countsHE) / 10000

#Deconvoluted
szfCluster <- quickCluster(countsHE)
szf_alClust <- normalizeBySums(countsHE, cluster=szfCluster, positive = TRUE)

# ---- SF-Difference ----
cellCols <- rainbow(nlevels(szfCluster))[as.integer(szfCluster)]
scaled_factors<- data.frame(DESeq=szf_HE/median(szf_HE), TMM=tmm/median(tmm), LibSize=libfactor/median(libfactor), Deconvolution=szf_alClust/median(szf_alClust), check.names=FALSE)

cairo_pdf("Macosko_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors, log="y", ylab="Size factors", cex.axis=1.5, cex.lab=1.8, xaxt='n')
text(labels=colnames(scaled_factors), xpd=TRUE, cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.1, srt=45, pos=2)
dev.off()

cairo_pdf("Macosko_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq, scaled_factors$Deconvolution, cex=0.6, pch=19, col="#00000073", xlab="DESeq", ylab="Deconvolution", xlim=c(0.1, 3.5), ylim=c(0.1, 20), log="xy", cex.axis=1.5, cex.lab=1.8)
abline(0, 1, col="dodgerblue")
dev.off()

cairo_pdf("Macosko_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$LibSize, scaled_factors$Deconvolution, cex=0.6, pch=19, col=cellCols, xlab="Library size", ylab="Deconvolution", xlim=c(0.1, 20), ylim=c(0.1, 20), log="xy", cex.axis=1.5, cex.lab=1.8)
abline(0, 1, col="dodgerblue")
dev.off()

cairo_pdf("Macosko_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM, scaled_factors$Deconvolution, cex=0.6, pch=19, col="#00000073", xlab="TMM", ylab="Deconvolution", xlim=c(0.1, 20), ylim=c(0.1, 20), log="xy", cex.axis=1.5, cex.lab=1.8)
abline(0, 1, col="dodgerblue")
dev.off()

distM <- as.dist(1-cor(countsHE, method="spearman"))
tsn <- tsne(distM)
plot(tsn,col=cellCols)

counts.pca <- prcomp(t(log(countsHE+1)))
plot(counts.pca$x[,1:2], col=cellCols, pch=19, cex=0.6)
tsn <- tsne(counts.pca$x[,1:10], initial_dims=10, k=2)
cairo_pdf("TSNE-Plot.pdf")
plot(tsn,col=cellCols, pch=19, cex=0.6)
dev.off()
