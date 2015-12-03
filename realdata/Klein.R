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
esd2 <- read.csv("../data/GSM1599497_ES_d2_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  
esd4 <- read.csv("../data/GSM1599498_ES_d4_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  
esd7 <- read.csv("../data/GSM1599499_ES_d7_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  

colnames(esd0) <- c("feature",c(1:933))
colnames(esd2) <- c("feature",c(934:1236))
colnames(esd4) <- c("feature",c(1237:1919))
colnames(esd7) <- c("feature",c(1920:2717))

counts <- merge (esd0,esd2,by="feature",all=T)
counts <- merge (counts, esd4, by="feature", all=T)
counts <- merge (counts, esd7, by="feature", all=T)
rownames(counts) <- counts$feature
counts <- counts[,!(colnames(counts) %in% "feature")]

celltype<- c(rep("d0",933),rep("d2", 303),rep("d4", 683),rep("d7" , 798))

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

#SF For AK method
szf_al <- normalizeBySums(countsHE)

#SF for AK with Cluster
szfCluster <- quickCluster(countsHE)
szf_alClust <- normalizeBySums(countsHE, clusters = szfCluster)


# ---- SF-Difference ----

all_sizefactors <- data.frame("SF" = szf_HE,"SF_after_Sum." = szf_al, "SF_after_SumClust" = szf_alClust)
boxplot(all_sizefactors,main="Distribution of all SF")
plot(szf_HE,szf_al,main="SF vs SF After Summation",xlab="Normal SF",ylab="SF after Summation",cex=0.6,pch=19,col="#00000073")
plot(szf_HE,szf_alClust,main="SF vs SF After Summation and Clustering",xlab="Normal SF",ylab="SF after Summation and Clustering",cex=0.6,pch=19,col="#00000073")
plot(szf_al,szf_alClust,main="SF After Summation vs SF After Summation and Clustering",xlab="SF after Summation",ylab="SF after Summation and Clustering",cex=0.6,pch=19,col="#00000073")

# ---- Normalize ----
countsHENorm <- as.data.frame(t(t(countsHE) / szf_HE))
countsHENorm_al <- as.data.frame(t(t(countsHE) / szf_al))

featuresNorm <- featureCalc(countsHENorm, htseq = FALSE, exprmin = 1)
featuresNorm_al <- featureCalc(countsHENorm_al, htseq = FALSE, exprmin = 1)
