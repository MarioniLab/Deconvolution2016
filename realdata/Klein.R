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

# ---- Normalize ----

countsHENorm <- as.data.frame(t(t(countsHE) / szf_HE))
countsHENorm_alClust <- as.data.frame(t(t(countsHE) / szf_alClust))


# ---- SF-Difference ----
scaled_factors<- data.frame("DESeq" = szf_HE / median(szf_HE),"TMM" = tmm/median(tmm), "Library size" = libfactor/median(libfactor), "Deconvolution" = szf_alClust/median(szf_alClust),check.names=FALSE)
boxplot(scaled_factors,log="y",ylab="Scaled Normalizationfactor",cex.axis=1.5,cex.lab=1.8)
plot(scaled_factors$Sizefactor,scaled_factors$Deconvoluted,cex=0.6,pch=19,col="#00000073",xlab="DESeq",ylab="Deconvolution",xlim=c(0.1,3.5),ylim=c(0.1,3.5),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")

plot(scaled_factors$LibrarySize,scaled_factors$Deconvoluted,cex=0.6,pch=19,col="#00000073",xlab="Librarysize",ylab="Deconvolution",xlim=c(0.1,8),ylim=c(0.1,8),log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")

# ---- Diferential-Expression ----

## Subset Data for EdgeR
cells <- data.frame("Cells" = celltype)

# Setup EdgeR
desgn <- model.matrix(~factor(cells$Cells))
y <- DGEList(countsHE)

#Sizefactors
y.sf <- y
y.sf$samples$norm.factors <- szf_HE/y.sf$samples$lib.size
y.sf <- estimateDisp(y.sf,desgn)
fit.sf <- glmFit(y.sf,desgn)
res.sf <- glmLRT(fit.sf)

#TMM
y.tmm <- y
y.tmm$samples$norm.factors <- tmm/y.tmm$samples$lib.size
y.tmm <- estimateDisp(y.tmm,desgn)
fit.tmm <- glmFit(y.tmm,desgn)
res.tmm <- glmLRT(fit.tmm)

#Library Size
y.lib <- y
y.lib$samples$norm.factors <- 1
y.lib <- estimateDisp(y.lib,desgn)
fit.lib <- glmFit(y.lib,desgn)
res.lib <- glmLRT(fit.lib)

#Deconvoluted
y.al <- y
y.al$samples$norm.factors <- szf_alClust/y.al$samples$lib.size
y.al <- estimateDisp(y.al,desgn)
fit.al <- glmFit(y.al,desgn)
res.al <- glmLRT(fit.al)

## Comparison

x.sf <- decideTestsDGE(res.sf)
x.lib <- decideTestsDGE(res.lib)
x.tmm <- decideTestsDGE(res.tmm)
x.al <- decideTestsDGE(res.al)


out.file <- "Klein_output.txt"
write.table(file=out.file, data.frame(Method="Size factor", Total=sum(x.sf!=0), Down=sum(x.sf<0), Up=sum(x.sf>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(file=out.file, data.frame(Method="TMM", Total=sum(x.tmm!=0), Down=sum(x.tmm<0), Up=sum(x.tmm>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="Library size", Total=sum(x.lib!=0), Down=sum(x.lib<0), Up=sum(x.lib>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="Deconvolution", Total=sum(x.al!=0), Down=sum(x.al<0), Up=sum(x.al>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

write.table(file=out.file, data.frame(Method="vs SF", Total=sum(x.sf!=0 & x.al!=0), Down=sum(x.sf<0 & x.al < 0), Up=sum(x.sf>0 & x.al > 0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="vs TMM", Total=sum(x.tmm!=0 & x.al!=0), Down=sum(x.tmm<0 & x.al < 0), Up=sum(x.tmm>0 & x.al > 0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="vs lib", Total=sum(x.lib!=0 & x.al!=0), Down=sum(x.lib<0 & x.al < 0), Up=sum(x.lib>0 & x.al > 0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
