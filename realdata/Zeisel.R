# AK47 Normalization

# ---- Packages ----
pkgs <- "/nfs/research2/marioni/Karsten/utils/Rpack/"
library('rJava',lib.loc=pkgs)
library('xlsxjars',lib.loc=pkgs)
library('xlsx',lib.loc=pkgs)
library('DESeq2')
library('plyr')
library('edgeR')
library('dynamicTreeCut',lib.loc=pkgs)
library('moduleColor',lib.loc=pkgs)
library('lattice',lib.loc=pkgs)
library('ggplot2',lib.loc=pkgs)
library('gridExtra',lib.loc=pkgs)
library('AK47')
source("./functions.R") 
source("./DM.R")

# ---- Data ----
counts <- read.delim("../data/expression_mRNA_17-Aug-2014.txt", stringsAsFactors = FALSE, header =FALSE)
molecules <- read.delim("../data/SilverBulletCTRLConc.txt",stringsAsFactors = FALSE,header = TRUE)
spikecounts <- read.delim("../data/expression_spikes_17-Aug-2014.txt",stringsAsFactors = FALSE, header=FALSE)
mitcounts <- read.delim("../data/expression_mito_17-Aug-2014.txt", stringsAsFactors = FALSE, header= FALSE)
# Organize Data
rownames(counts)[12:nrow(counts)] <- counts$V1[12:nrow(counts)] # Genenames as rownames 
rownames(counts)[1:10] <- counts$V2[1:10] # Cell metadata type as rownames
counts <- counts[,c(-1,-2)]            # remove the coloumns that are now rownames
colnames(counts) <- as.character(counts["cell_id",]) # make cell_id colname
cells <- data.frame(t(counts[1:10,]))  # extract cell metdata into seperate dataframe
colnames(cells)[colnames(cells)=="group.."] <- "group"
counts <- counts[12:nrow(counts),]     # remove metadata from countdata
counts <- as.matrix(counts)            # change typeof coloumns to numeric 
mode(counts) <- "numeric"
counts <- as.data.frame(counts)

# Organize spike data (analogous to count data)
rownames(spikecounts)[11:nrow(spikecounts)] <- spikecounts$V1[11:nrow(spikecounts)]
rownames(spikecounts)[1:10] <- spikecounts$V2[1:10]
spikecounts <- spikecounts[,c(-1,-2)]
colnames(spikecounts) <- spikecounts["cell_id",]
spikecounts <- spikecounts[12:nrow(spikecounts),]
spikecounts <- as.matrix(spikecounts)
mode(spikecounts) <- "numeric"
spikecounts <- as.data.frame(spikecounts)

# Organize mito data (analogous to count data)
rownames(mitcounts)[11:nrow(mitcounts)] <- mitcounts$V1[11:nrow(mitcounts)]
rownames(mitcounts)[1:10] <- mitcounts$V2[1:10]
mitcounts <- mitcounts[,c(-1,-2)]
colnames(mitcounts) <- mitcounts["cell_id",]
mitcounts <- mitcounts[12:nrow(mitcounts),]
mitcounts <- as.matrix(mitcounts)
mode(mitcounts) <- "numeric"
mitcounts <- as.data.frame(mitcounts)
mitcounts <- mitcounts[,colnames(counts)]

# Combine data sets into one expression matrix 
counts <- rbind(counts,mitcounts,spikecounts)


minExpr <- 1

features <- featureCalc(counts = counts,htseq = FALSE, exprmin = minExpr) # calcuate feature data 
features$mito <- grepl("mt-",rownames(counts)) # Add logical coloumn with mt gene indices

# ---- Gene-Filter ----
# Remove lowly expressed Genes
IncCrit <- ncol(counts)/5

highE<- rowSums(counts) >= IncCrit
countsHE <- counts[highE,]

featuresHE <- featureCalc(countsHE, htseq = FALSE, minExpr)

# ---- Size-Factors ----
countsHEnoSpike <- countsHE[!featuresHE$spike,]

#Normal SF
geoMeans_HE <- apply(countsHEnoSpike, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))

szf_HE <- estimateSizeFactorsForMatrix(countsHEnoSpike, geoMeans = geoMeans_HE)

#Normal TMM
tmm <- calcNormFactors(as.matrix(countsHEnoSpike),method="TMM") * colSums(countsHEnoSpike)

#Normal Library Size
libfactor <- colSums(countsHEnoSpike) / 10000

#ERCC Size factors
erccfactor <- colSums(countsHE[featuresHE$spike,]) / 1000

#Deconvoluted
szfCluster <- quickCluster(countsHEnoSpike)
szf_alClust <- normalizeBySums(countsHEnoSpike, clusters = szfCluster)

# ---- SF-Difference ----
set.seed(1)
shffl <- sample(ncol(countsHE))
scaled_factors <- data.frame("DESeq" = (szf_HE/median(szf_HE))[shffl],"TMM" = (tmm/median(tmm))[shffl], "Library size" = (libfactor/median(libfactor))[shffl], "Deconvolution" = (szf_alClust/median(szf_alClust))[shffl],check.names=FALSE)
#Comparison of SF Distribution
#Colors
line.col <- "red"
olig.col <- rgb(1,0.5,0)
pyr.col <- "dodgerblue"
cell.col <- rgb(0,0,0) 

setCols <- sapply(cells$level1class[shffl], function(x) if (x == "oligodendrocytes") olig.col  else if (x == "pyramidal CA1") pyr.col else cell.col)

cairo_pdf("Zeisel_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors,log="y",ylab="Size factors",cex.axis=1.5,cex.lab=1.8,xaxt='n')
text(labels=colnames(scaled_factors),xpd=TRUE,cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.1,srt=45, pos=2)
dev.off()
#SF vs Deconvoluted

cairo_pdf("Zeisel_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,scaled_factors$Deconvolution,xlim=c(0.1,6),ylim=c(0.1,6),pch=16,xlab="DESeq",ylab="Deconvolution",log="xy",cex.axis=1.5,cex.lab=1.8,col=setCols)
abline(0,1,col=line.col,lwd=2)
dev.off()

#Libfactor vs Deconvoluted
cairo_pdf("Zeisel_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors[,"Library size"],scaled_factors$Deconvolution,xlim=c(0.1,7),ylim=c(0.1,7),pch=16,xlab="Library size",ylab="Deconvolution",log="xy",cex.axis=1.5,cex.lab=1.8,col=setCols)
abline(0,1,col=line.col,lwd=2)
dev.off()

#TMM vs Deconvoluted
cairo_pdf("Zeisel_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM,scaled_factors$Deconvolution,xlim=c(0.1,7),ylim=c(0.1,7),pch=16,col=setCols,xlab="TMM",ylab="Deconvolution",cex.axis=1.5,cex.lab=1.8,log="xy")
abline(0,1,col=line.col,lwd=2)
dev.off()

# LFC in sizefactors between sziefactor and deconvolution for the a) all b) pyramidal c) oligodendrocytes
# scaled_lib <- libfactor/median(libfactor[cells$group != "4"])
# scaled_szf <- szf_alClust/median(szf_alClust[cells$group != "4"])
# 
# lfc_al <- log2(scaled_lib/scaled_szf)
# lfc_oligo <- log2(scaled_lib[cells$group == "4"] / scaled_szf[cells$group == "4"])
# lfc_pyramidal <- log2(scaled_lib[cells$group == "3"] / scaled_szf[cells$group == "3"])
# lfc_micro <- log2(scaled_lib[cells$group == "5"] / scaled_szf[cells$group == "5"])
# 
# boxplot(lfc_al,lfc_oligo,lfc_pyramidal,lfc_micro)

# ---- Spike-Ins ----
cairo_pdf("Zeisel_ERCCvSF.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,erccfactor/median(erccfactor),ylab="Spike-in",xlab="DESeq",log="x",pch=16,col="#00000073",cex.axis=1.5,cex.lab=1.8)
dev.off()
cairo_pdf("Zeisel_ERCCvTMM.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM,erccfactor/median(erccfactor),ylab="Spike-in",xlab="TMM",log="x",cex=0.6,pch=19,col="#00000073",cex.axis=1.5,cex.lab=1.8)
dev.off()
cairo_pdf("Zeisel_ERCCvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$Deconvolution,erccfactor/median(erccfactor),ylab="Spike-in",xlab="Deconvolution",log="x",cex=0.6,pch=19,col="#00000073",cex.axis=1.5,cex.lab=1.8)
dev.off()
cairo_pdf("Zeisel_ERCCvLib.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors[,"Library size"],erccfactor/median(erccfactor),ylab="Spike-in",xlab="Library size",log="x",cex=0.6,pch=19,col="#00000073",cex.axis=1.5,cex.lab=1.8)
dev.off()

# ---- Normalization ----

counts.sf <- as.data.frame(t(t(countsHEnoSpike) / szf_HE))
counts.tmm <- as.data.frame(t(t(countsHEnoSpike) / tmm))
counts.lib <- as.data.frame(t(t(countsHEnoSpike) / libfactor))
counts.decon <- as.data.frame(t(t(countsHEnoSpike) / szf_alClust))

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
    write.table(comparisonMatrix,paste("HVGcomparisonZeisel",x,".txt",sep=""))
}
