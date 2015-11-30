## Clustering a binarized set of lowly expressed genes achieves similar accuracy while allowing better DE calling


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
library('dynamicTreeCut',lib.loc=pkgs)
library('moduleColor',lib.loc=pkgs)
library('lattice',lib.loc=pkgs)
library('ggplot2',lib.loc=pkgs)
library('gridExtra',lib.loc=pkgs)
source("./functions.R") 
source("./DM.R")
source("../prenormSum/sum_code.R")
# ---- Data ----
counts <- read.delim("../data/expression_mRNA_17-Aug-2014.txt", stringsAsFactors = FALSE, header =FALSE)
molecules <- read.delim("../data/SilverBulletCTRLConc.txt",stringsAsFactors = FALSE,header = FALSE)
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

# Create a neuronset and a normal set
# Not used at the moment
neurons <- c("^CA1","^CA2","^Int","^S1")
nset <- grepl(paste(neurons,collapse="|"),cells$level2class)
ncounts <- counts[,nset]
ncells <- cells[nset,]
ncells <- droplevels(ncells)

minExpr <- 1

features <- featureCalc(counts = counts,htseq = FALSE, exprmin = minExpr) # calcuate feature data 
features$mito <- grepl("mt-",rownames(counts)) # Add logical coloumn with mt gene indices

nfeatures<- featureCalc(counts= ncounts,htseq = FALSE, exprmin = minExpr) # calcuate feature data 
nfeatures$mito <- grepl("mt-",rownames(ncounts)) # Add logical coloumn with mt gene indices

# ---- Gene-Filter ----
# Dataset is split into HE and LE, LE is everything below a mean of 0.2 
# Create data set for HEG

IncCrit <- ncol(counts)/5

highE<- rowSums(counts) >= IncCrit
countsHE <- counts[highE,]

featuresHE <- featureCalc(countsHE, htseq = FALSE, minExpr)

# Create data set for LEG
lowE <- rowSums(counts) < IncCrit & rowSums(counts) > 0

countsLE <- counts[lowE,]


#Binarize dataset for LEG

countsLE <- data.frame((countsLE > 0) * 1)

# ---- LEG-Cluster ----

dismLE <- dist(t(countsLE),method='binary')
clustLE <- hclust(dismLE,method = "ward.D2")


cutTLE <- cutreeDynamic(clustLE, minClusterSize = 30, method = 'hybrid',
                      distM =as.matrix(dismLE), deepSplit = 1, pamStage = TRUE, respectSmallClusters = TRUE, verbose = 2)

colLE <- labels2colors(cutTLE)


cutTLE_deep <- cutreeDynamic(clustLE, minClusterSize = 10, method = 'hybrid',
                      distM =as.matrix(dismLE), deepSplit = 4, pamStage = TRUE, respectSmallClusters = TRUE, verbose = 2)

colLE_deep <- labels2colors(cutTLE_deep)

kmns <- kmeans(dismLE,centers = 300)


# par(mfrow=c(1,2))
# plot(clustLE,labels = FALSE,main="Dendogram LE ")
# plotgroups(clustLE,colLE,as.numeric(cells$group))
# plot(clustLE,labels = FALSE,main="Dendogram LE ")
# plotgroups(clustLE,colLE_deep,as.numeric(cells$level2class))

# ---- Subgroup ----

set.seed(1)
sumgroups <- ClusterSums(countsHE,factor(cutTLE), 10)
sumgroups_deep <- ClusterSums(countsHE,factor(cutTLE_deep),10)
sumgroups_kmns <- ClusterSums(countsHE,factor(unname(kmns$cluster)),max(kmns$size))


# ---- Size-Factors ----
# SF for the summed count matrices
geoMeans <- apply(sumgroups[[1]], 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
sf <- estimateSizeFactorsForMatrix(sumgroups[[1]],geoMeans = geoMeans)

geoMeans_deep <- apply(sumgroups_deep[[1]], 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
sf_deep <- estimateSizeFactorsForMatrix(sumgroups_deep[[1]],geoMeans = geoMeans_deep)

geoMeans_kmns <- apply(sumgroups_kmns[[1]], 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
sf_kmns <- estimateSizeFactorsForMatrix(sumgroups_kmns[[1]],geoMeans = geoMeans_kmns)

# Create sizefactor vector based on the summed count matrix for countHE
szf <- unlist(sapply(names(sf),function(sumgrp) rep(sf[sumgrp],length(sumgroups[[2]][[sumgrp]]))))
names(szf)<-unname(unlist(sapply(names(sf), function(sumgrp) sumgroups[[2]][[sumgrp]])))

szf_deep <- unlist(sapply(names(sf_deep),function(sumgrp) rep(sf_deep[sumgrp],length(sumgroups_deep[[2]][[sumgrp]]))))
names(szf_deep)<-unname(unlist(sapply(names(sf_deep), function(sumgrp) sumgroups_deep[[2]][[sumgrp]])))

szf_kmns <- unlist(sapply(names(sf_kmns),function(sumgrp) rep(sf_kmns[sumgrp],length(sumgroups_kmns[[2]][[sumgrp]]))))
names(szf_kmns)<-unname(unlist(sapply(names(sf_kmns), function(sumgrp) sumgroups_kmns[[2]][[sumgrp]])))

# SF for countsHE 
geoMeans_HE <- apply(countsHE, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
szf_HE <- estimateSizeFactorsForMatrix(countsHE,geoMeans = geoMeans_HE)

#Sf For Aarons method
szf_al <- normalizeBySums(countsHE)

# ---- SF-Difference ----
sfdiff <- sapply(names(sf), function(x) szf_HE[sumgroups[[2]][[x]]] - sf[x])
sfdiff_deep <- sapply(names(sf_deep), function(x) szf_HE[sumgroups_deep[[2]][[x]]] - sf_deep[x])
sfdiff_kmns <- sapply(names(sf_kmns), function(x) szf_HE[sumgroups_kmns[[2]][[x]]] - sf_kmns[x])

boxplot(sfdiff,main="Summing after Hclust (k=10)")
boxplot(sfdiff_deep,main="Summing after Hclust (k=58)")
boxplot(sfdiff_kmns,main="Summing after Kmeans (k=300)")

all_sizefactors <- data.frame(szf_HE,szf_al,szf_kmns)
boxplot(all_sizefactors,main="Distribution of all SF")
pairs(all_sizefactors,cex=0.4,main="Pairwise Comparison of SF")

# ---- Normalize ----
countsHENorm <- as.data.frame(t(t(countsHE) / szf_HE))
countsHENorm_sum <- as.data.frame(t(t(countsHE) /szf[colnames(countsHE)]))
countsHENorm_sum_deep <- as.data.frame(t(t(countsHE) /szf_deep[colnames(countsHE)]))
countsHENorm_sum_kmns <- as.data.frame(t(t(countsHE) /szf_kmns[colnames(countsHE)]))
countsHENorm_al <- as.data.frame(t(t(countsHE) / szf_al))

featuresNorm <- featureCalc(countsHENorm, htseq = FALSE, exprmin = 1)
featuresNorm_sum <- featureCalc(countsHENorm_sum, htseq = FALSE, exprmin = 1)
featuresNorm_sum_deep <- featureCalc(countsHENorm_sum_deep, htseq = FALSE, exprmin = 1)
featuresNorm_sum_kmns <- featureCalc(countsHENorm_sum_kmns, htseq = FALSE, exprmin = 1)
featuresNorm_al <- featureCalc(countsHENorm_al, htseq = FALSE, exprmin = 1)

# ---- Evaluate-Noise ----
par(mfrow=c(3,2))
spike_colNorm <- as.numeric(featuresNorm$spike) + 1
plot(featuresNorm$mean,featuresNorm$cv2,log="xy",cex=0.3,col=spike_colNorm,pch=19,main="Not Summed")
curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")


spike_colNorm_sum <- as.numeric(featuresNorm_sum$spike) + 1
plot(featuresNorm_sum$mean,featuresNorm_sum$cv2,log="xy",cex=0.3,col=spike_colNorm_sum,pch=19,main="Summed Hclust(k=10)")
curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")


spike_colNorm_sum_deep <- as.numeric(featuresNorm_sum_deep$spike) + 1
plot(featuresNorm_sum_deep$mean,featuresNorm_sum_deep$cv2,log="xy",cex=0.3,col=spike_colNorm_sum_deep,pch=19,main="Summed Hclust(k=58)")
curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")


spike_colNorm_sum_kmns <- as.numeric(featuresNorm_sum_kmns$spike) + 1
plot(featuresNorm_sum_kmns$mean,featuresNorm_sum_kmns$cv2,log="xy",cex=0.3,col=spike_colNorm_sum_kmns,pch=19,main="Summed Kmeans(k=300)")
curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")


spike_colNorm_al <- as.numeric(featuresNorm_al$spike) + 1
plot(featuresNorm_al$mean,featuresNorm_al$cv2,log="xy",cex=0.3,col=spike_colNorm_al,pch=19,main="Summed method")
curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")
dev.off()

dm <- DM(meanGenes=featuresHE$mean,CV2Genes=featuresHE$cv2)
dmNorm <- DM(meanGenes=featuresNorm$mean,CV2Genes=featuresNorm$cv2)
dmNorm_sum <- DM(meanGenes=featuresNorm_sum$mean,CV2Genes=featuresNorm_sum$cv2)
dmNorm_sum_deep <- DM(meanGenes=featuresNorm_sum_deep$mean,CV2Genes=featuresNorm_sum_deep$cv2)
dmNorm_sum_kmns <- DM(meanGenes=featuresNorm_sum_kmns$mean,CV2Genes=featuresNorm_sum_kmns$cv2)
dmNorm_al <- DM(meanGenes=featuresNorm_al$mean,CV2Genes=featuresNorm_al$cv2)

dmERCC <- dm[grepl("ERCC",rownames(countsHE))]
dmERCC_norm <- dmNorm[grepl("ERCC",rownames(countsHE))]
dmERCC_sum <- dmNorm_sum[grepl("ERCC",rownames(countsHE))]
dmERCC_sum_deep <- dmNorm_sum_deep[grepl("ERCC",rownames(countsHE))]
dmERCC_sum_kmns <- dmNorm_sum_kmns[grepl("ERCC",rownames(countsHE))]
dmERCC_al <- dmNorm_al[grepl("ERCC",rownames(countsHE))]

boxplot(data.frame(dmERCC,dmERCC_norm,dmERCC_sum,dmERCC_sum_deep,dmERCC_sum_kmns,dmERCC_al),main="Boxplot of DM of ERCCs")
# ---- SystematicClusterDifference ----

## First Plot a Boxplot of the numbers of Zeros per Cluster
nzeros <- data.frame(t(apply(countsHE,2, function(x) sum(x == 0))))
rownames(nzeros) <- "NZeros"
plotcutTLE<- drawBoxplot(nzeros,cutTLE,"NZeros")
plotcutTLE_deep<- drawBoxplot(nzeros,cutTLE_deep,"NZeros")
plot_kmns<- drawBoxplot(nzeros,kmns$cluster,"NZeros")

plot(plotcutTLE)
plot(plotcutTLE_deep)
plot(plot_kmns)

#Plot the numbers of zeros per Cluster to the difference of the per Cluster median of the per Sum median difference to szf_HE

nzeros_clust <- data.frame(t(nzeros),"Cluster" = factor(cutTLE))
nzeros_deep <- data.frame(t(nzeros),"Cluster" = factor(cutTLE_deep))
nzeros_kmns <- data.frame(t(nzeros),"Cluster" = factor(kmns$cluster))

##First Calculate meidan of Zeros per Cluster
nzeros_clust_medians <- sapply(levels(nzeros_clust$Cluster), function(x) median(nzeros_clust[nzeros_clust$Cluster == x,"NZeros"]))
nzeros_deep_medians <- sapply(levels(nzeros_deep$Cluster), function(x) median(nzeros_deep[nzeros_deep$Cluster == x,"NZeros"]))
nzeros_kmns_medians <- sapply(levels(nzeros_kmns$Cluster), function(x) median(nzeros_kmns[nzeros_kmns$Cluster == x,"NZeros"]))


## Second Calculate per Sum median of differences to szf_HE
perSumDiff_clust <- unlist(lapply(sfdiff,median))
perSumDiff_deep <- unlist(lapply(sfdiff_deep,median))
perSumDiff_kmns <- unlist(lapply(sfdiff_kmns,median))

#Remove everything but the cluster name from the name
names(perSumDiff_clust) <- sapply(strsplit(names(perSumDiff_clust), "_"), "[",2)
names(perSumDiff_deep) <- sapply(strsplit(names(perSumDiff_deep), "_"), "[",2)
names(perSumDiff_kmns) <- sapply(strsplit(names(perSumDiff_kmns), "_"), "[",2)
                               
## Third calc the per cluster median of the per sum median of the differences to szf_HE
perClustDiff_clust <- sapply(levels(factor(as.numeric(names(perSumDiff_clust)))), function(x) median(perSumDiff_clust[names(perSumDiff_clust) == x]))
perClustDiff_deep <- sapply(levels(factor(as.numeric(names(perSumDiff_deep)))), function(x) median(perSumDiff_deep[names(perSumDiff_deep) == x]))
perClustDiff_kmns <- sapply(levels(factor(as.numeric(names(perSumDiff_kmns)))), function(x) median(perSumDiff_kmns[names(perSumDiff_kmns) == x]))

# Plot
par(mfrow=c(2,2),pch=19,cex=0.5,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
plot(nzeros_clust_medians,perClustDiff_clust,xlab="Median #Zeros per Cluster",ylab="PerClust Median of median per Sum of SFClust - SF",main="HClust (k=10)")
plot(nzeros_deep_medians,perClustDiff_deep,xlab="Median #Zeros per Cluster",ylab="PerClust Median of median per Sum of SFClust - SF",main="HClust (k=58)")
plot(nzeros_kmns_medians,perClustDiff_kmns,xlab="Median #Zeros per Cluster",ylab="PerClust Median of median per Sum of SFClust - SF",main="Kmeans (k=300)")
plot(t(nzeros)[,1],szf_HE,xlab="Cell-#Zeros",ylab="Cell-Sizefactors",main="Without Summation")

# 
# ---- ExtremeSFValues ----
par(mfrow=c(1,2))
plot(szf_al,t(nzeros)[,1]/median(t(nzeros)[,1]),xlab="Sizefactors After Summation",ylab="Relative Amount of Zeros")
plot(szf_HE,t(nzeros)[,1]/median(t(nzeros)[,1]),xlab="Sizefactors",ylab="Relative Amount of Zeros")

par(mfrow=c(1,2))
plot(szf_al,colSums(countsHE),xlab="Sizefactors After Summation",ylab="Relative Lib.Size")
plot(szf_HE,colSums(countsHE),xlab="Sizefactors",ylab="Relative Lib.Size")

# ---- Spike-Ins ----
par(mfrow=c(1,2))
plot(szf_HE,colSums(countsHE[!featuresHE$spike,])/colSums(countsHE[featuresHE$spike,]),ylab="Cellsize as Cellreads / ERCCreads",xlab="Sizefactors",log="x")
plot(szf_al,colSums(countsHE[!featuresHE$spike,])/colSums(countsHE[featuresHE$spike,]),ylab="Cellsize as Cellreads / ERCCreads",xlab="Sizefactors after Summation",log="x")

par(mfrow=c(1,2))
plot(szf_HE,colSums(countsHE[featuresHE$spike,]),ylab=" ERCCreads",xlab="Sizefactors",log="x")
plot(szf_al,colSums(countsHE[featuresHE$spike,]),ylab=" ERCCreads",xlab="Sizefactors after Summation",log="x")

par(mfrow=c(1,2))
plot(colSums(countsHE[featuresHE$spike,]),colSums(countsHE[!featuresHE$spike,])/colSums(countsHE[featuresHE$spike,]),ylab="Cellsize as Cellreads / ERCCreads",xlab="ERCC reads",log="x")
plot(colSums(countsHE[featuresHE$spike,]),colSums(countsHE[!featuresHE$spike,])/colSums(countsHE[featuresHE$spike,]),ylab="Cellsize as Cellreads / ERCCreads",xlab="ERCC reads",log="x")
# ---- Effect-on-HVG ---

##Order of Variability
plot(order(dm),order(dmNorm),col=grepl("ERCC",rownames(countsHE)) + 1,cex=0.5,pch=19)
plot(order(dm),order(dmNorm_al),col=grepl("ERCC",rownames(countsHE)) + 1,cex=0.5,pch=19)
plot(order(dmNorm),order(dmNorm_al),col=grepl("ERCC",rownames(countsHE)) + 1,cex=0.5,pch=19)

featuresHE_Ord<- featuresHE[order(-dm),]
featuresNorm_Ord<- featuresNorm[order(-dmNorm),]
featuresNorm_alOrd <- featuresNorm_al[order(-dmNorm_al),]

featuresHE$hvg <- rownames(featuresHE) %in% rownames(featuresHE_Ord)[1:1000]
featuresNorm$hvg <- rownames(featuresNorm) %in% rownames(featuresNorm_Ord)[1:1000]
featuresNorm_al$hvg <- rownames(featuresNorm_al) %in% rownames(featuresNorm_alOrd)[1:1000]

venHE <- which(featuresHE$hvg)
venNorm <- which(featuresNorm$hvg)
venNorm_al <- which(featuresNorm_al$hvg)

venn.diagram(list("Non-Norm" = venHE,"Sizefacotrs" = venNorm, "Summation" = venNorm_al), fill = c("red","green","blue"),
             alpha= c(0.5,0.5,0.5), cex = 1 , cat.fontface = 4, fontfamily = 3,imagetype="png",
             filename = "VenDiagram2.png")

# ---- Diferential-Expression ----

sbst <- countsHE[,grepl("microglia|oligodendrocytes",cells$level1class)]
cellssbst <- cells[grepl("microglia|oligodendrocytes",cells$level1class),]
cellssbst$sizeFactor <- szf_HE[grepl("microglia|oligodendrocytes",cells$level1class)]

cellssbst_al <- cells[grepl("microglia|oligodendrocytes",cells$level1class),]
cellssbst_al$sizeFactor <- szf_al[grepl("microglia|oligodendrocytes",cells$level1class)]


## Load into DESeq2
dds_sbst <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst, design =~ level1class)
dds_sbst <- estimateDispersions(dds_sbst)
dds_sbst <- nbinomWaldTest(dds_sbst)


dds_sbst_al <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst_al, design =~ level1class)
dds_sbst_al <- estimateDispersions(dds_sbst_al)
dds_sbst_al <- nbinomWaldTest(dds_sbst_al)
