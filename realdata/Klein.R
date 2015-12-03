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
scaled_factors<- data.frame("Sizefactor" = szf_HE / median(szf_HE),"TMM" = tmm/median(tmm), "LibrarySize" = libfactor/median(libfactor), "Deconvoluted" = szf_alClust/median(szf_alClust))
boxplot(scaled_factors,log="y",ylab="Scaled Normalizationfactor")
plot(scaled_factors$Sizefactor,scaled_factors$Deconvoluted,cex=0.6,pch=19,col="#00000073",xlab="Sizefactor",ylab="Deconvoluted",xlim=c(0.1,3.5),ylim=c(0.1,3.5),log="xy")
abline(0,1,col="dodgerblue")

plot(scaled_factors$LibrarySize,scaled_factors$Deconvoluted,cex=0.6,pch=19,col="#00000073",xlab="Librarysize",ylab="Deconvoluted",xlim=c(0.1,8),ylim=c(0.1,8),log="xy")
abline(0,1,col="dodgerblue")

# ---- Diferential-Expression ----

#Construct colData for both Sets
colData_szf <- data.frame("Cells" = celltype, "sizeFactor" = szf_HE)
colData_szf_alClust <- data.frame("Cells" = celltype, "sizeFactor" = szf_alClust)
colData_lib <- data.frame("Cells" = celltype, "sizeFactor" = libfactor)
colData_tmm <- data.frame("Cells" = celltype, "sizeFactor" = tmm)


## Load into DESeq2 ! Careful takes a lot of time !
dds_szf<- DESeqDataSetFromMatrix(countData = countsHE, colData = colData_szf, design = ~ Cells)
dds_szf <- estimateDispersions(dds_szf)
dds_szf <- nbinomWaldTest(dds_szf)
rslt_szf <- results(dds_szf)
write.csv(as.data.frame(rslt_szf),"Klein_Results_SF.csv")

dds_szf_al <- DESeqDataSetFromMatrix(countData = countsHE, colData =colData_szf_alClust, design = ~ Cells)
dds_szf_al <- estimateDispersions(dds_szf_al)
dds_szf_al <- nbinomWaldTest(dds_szf_al)
rslt_szf_al <- results(dds_szf_al)
write.csv(as.data.frame(rslt_szf_al),"Klein_Results_Decon.csv")

dds_lib <- DESeqDataSetFromMatrix(countData = countsHE, colData =colData_lib, design = ~ Cells)
dds_lib <- estimateDispersions(dds_lib)
dds_lib <- nbinomWaldTest(dds_lib)
rslt_lib <- results(dds_lib)
write.csv(as.data.frame(rslt_lib),"Klein_Results_Lib.csv")

dds_tmm <- DESeqDataSetFromMatrix(countData = countsHE, colData =colData_tmm, design = ~ Cells)
dds_tmm <- estimateDispersions(dds_tmm)
dds_tmm <- nbinomWaldTest(dds_tmm)
rslt_tmm <- results(dds_tmm)
write.csv(as.data.frame(rslt_tmm),"Klein_Results_TMM.csv")
