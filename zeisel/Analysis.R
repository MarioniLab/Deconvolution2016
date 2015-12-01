## AK47 Normalization

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
library('pheatmap',lib.loc=pkgs)
library("RColorBrewer")
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

# Create a neuronset and a normal set

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
szf_HE <- estimateSizeFactorsForMatrix(countsHEnoSpike,geoMeans = geoMeans_HE)

#SF For AK method
szf_al <- normalizeBySums(countsHEnoSpike)

#SF for AK with Cluster
szfCluster <- quickCluster(countsHEnoSpike)
szf_alClust <- normalizeBySums(countsHEnoSpike, clusters = szfCluster)


# ---- SF-Difference ----

all_sizefactors <- data.frame("SF" = szf_HE,"SF_after_Sum." = szf_al, "SF_after_SumClust" = szf_alClust)
boxplot(all_sizefactors,main="Distribution of all SF")
plot(szf_HE,szf_al,main="SF vs SF After Summation",xlab="Normal SF",ylab="SF after Summation",cex=0.6,pch=19,col="#00000073")
plot(szf_HE,szf_alClust,main="SF vs SF After Summation and Clustering",xlab="Normal SF",ylab="SF after Summation and Clustering",cex=0.6,pch=19,col="#00000073")
plot(szf_al,szf_alClust,main="SF After Summation vs SF After Summation and Clustering",xlab="SF after Summation",ylab="SF after Summation and Clustering",cex=0.6,pch=19,col="#00000073")

# ---- Normalize ----
countsHENorm <- as.data.frame(t(t(countsHEnoSpike) / szf_HE))
countsHENorm_al <- as.data.frame(t(t(countsHEnoSpike) / szf_al))

featuresNorm <- featureCalc(countsHENorm, htseq = FALSE, exprmin = 1)
featuresNorm_al <- featureCalc(countsHENorm_al, htseq = FALSE, exprmin = 1)

# ---- Evaluate-Noise ----
# spike_colNorm <- as.numeric(featuresNorm$spike) + 1
# plot(featuresNorm$mean,featuresNorm$cv2,log="xy",cex=0.3,col=spike_colNorm,pch=19,main="Not Summed")
# curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")
# 
# spike_colNorm_al <- as.numeric(featuresNorm_al$spike) + 1
# plot(featuresNorm_al$mean,featuresNorm_al$cv2,log="xy",cex=0.3,col=spike_colNorm_al,pch=19,main="Summed method")
# curve((sqrt(x)/x)^2,from=0.001,to=1000,add = T, col = "red")
# dev.off()
# 
# dm <- DM(meanGenes=featuresHE$mean,CV2Genes=featuresHE$cv2)
# dmNorm <- DM(meanGenes=featuresNorm$mean,CV2Genes=featuresNorm$cv2)
# dmNorm_al <- DM(meanGenes=featuresNorm_al$mean,CV2Genes=featuresNorm_al$cv2)
# 
# dmERCC <- dm[grepl("ERCC",rownames(countsHE))]
# dmERCC_norm <- dmNorm[grepl("ERCC",rownames(countsHE))]
# dmERCC_al <- dmNorm_al[grepl("ERCC",rownames(countsHE))]
# 
# boxplot(data.frame(dmERCC,dmERCC_norm,dmERCC_al),main="Boxplot of DM of ERCCs")

# ---- ExtremeSFValues ----

# par(mfrow=c(1,2))
# plot(szf_al,t(nzeros)[,1]/median(t(nzeros)[,1]),xlab="Sizefactors After Summation",ylab="Relative Amount of Zeros")
# plot(szf_HE,t(nzeros)[,1]/median(t(nzeros)[,1]),xlab="Sizefactors",ylab="Relative Amount of Zeros")

plot(szf_al,colSums(countsHE),xlab="Sizefactors After Summation",ylab="Total Lib.Size",cex=0.6,pch=19,col="#00000073")
plot(szf_HE,colSums(countsHE),xlab="Sizefactors",ylab="Total Lib.Size",cex=0.6,pch=19,col="#00000073")

# ---- Sizefactor-Information ----

plot(szf_HE,colSums(countsHE[!featuresHE$spike,])/colSums(countsHE[featuresHE$spike,]),ylab="Cellsize as Cellreads / ERCCreads",xlab="Sizefactors",cex=0.6,pch=19,col="#00000073",main="Cellsize vs. SF")
plot(szf_al,colSums(countsHE[!featuresHE$spike,])/colSums(countsHE[featuresHE$spike,]),ylab="Cellsize as Cellreads / ERCCreads",xlab="Sizefactors after Summation",cex=0.6,pch=19,col="#00000073",main="Cellsize vs SF after Sum")

plot(szf_HE,colSums(countsHE),ylab="Total Counts",xlab="Sizefactors",cex=0.6,pch=19,col="#00000073")
plot(szf_al,colSums(countsHE),ylab="Total Counts",xlab="Sizefactors after Summation",cex=0.6,pch=19,col="#00000073")

# ---- Spike-Ins ----

plot(szf_HE,colSums(countsHE[featuresHE$spike,]),ylab=" ERCCreads",xlab="Sizefactors",log="x",cex=0.6,pch=19,col="#00000073")
plot(szf_al,colSums(countsHE[featuresHE$spike,]),ylab=" ERCCreads",xlab="Sizefactors after Summation",log="x",cex=0.6,pch=19,col="#00000073")

## Compute conversionFactor
betafactor <- apply(countsHE[featuresHE$spike,],2, function(x) lm(x ~ molecules[molecules$ERCC_ID %in% rownames(featuresHE)[featuresHE$spike],"molecules_in_each_chamber"])$coefficients[2])

## Sizefactors dont correspond to conversion factors
plot(szf_HE,betafactor,cex=0.6,pch=19,col="#00000073",xlab="Sizefactors",ylab="Conversion Factors",main="Sizefactors are not assoicated with Conversion Factors")

# ---- Effect-on-HVG ---

##Order of Variability
plot(order(dm),order(dmNorm),cex=0.5,pch=19)
plot(order(dm),order(dmNorm_al),cex=0.5,pch=19)
plot(order(dmNorm),order(dmNorm_al),cex=0.5,pch=19)

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


## Load into DESeq2 ! Careful takes a lot of time !
dds_sbst <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst, design =~ level1class)
dds_sbst <- estimateDispersions(dds_sbst)
dds_sbst <- nbinomWaldTest(dds_sbst)
rslt_sbst <- results(dds_sbst)

dds_sbst_al <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst_al, design =~ level1class)
dds_sbst_al <- estimateDispersions(dds_sbst_al)
dds_sbst_al <- nbinomWaldTest(dds_sbst_al)

rslt_sbst_al <- results(dds_sbst_al)

up_sbst <- rownames(rslt_sbst)[rslt_sbst$log2FoldChange > 0 & rslt_sbst$padj < 0.1]
up_sbst_al <- rownames(rslt_sbst_al)[rslt_sbst_al$log2FoldChange > 0 & rslt_sbst_al$padj < 0.1]

down_sbst <- rownames(rslt_sbst)[rslt_sbst$log2FoldChange < 0 & rslt_sbst$padj < 0.1]
down_sbst_al <- rownames(rslt_sbst_al)[rslt_sbst_al$log2FoldChange < 0 & rslt_sbst$padj < 0.1]

int_up <- intersect(up_sbst,up_sbst_al)
int_down <- intersect(down_sbst, down_sbst_al)

int_up_down <- intersect(up_sbst,down_sbst_al)
int_down_up <- intersect(down_sbst,up_sbst_al)

barplot(c(length(up_sbst),length(up_sbst_al),length(int_up),length(down_sbst),length(down_sbst_al),length(int_down),length(int_up_down),length(int_down_up)),width=0.5,names.arg=c("SF","SF_AS","Overlap","SF","SF_AS","Overlap","Overlap","Overlap"),col=c("#FB6A4A","#FB6A4A","#EF3B2C","#6BAED6","#6BAED6","#2171B5","seagreen2","tomato2"),cex.names=0.7,space=c(0.1,0.1,0.1,0.5,0.1,0.1,0.5,0.1),main="Oligodendrocyte v. Microglia")
legend("topright",fill=c("#FB6A4A","#6BAED6","seagreen2","tomato2"),legend=c("Upregulated","Downregulated","Down after Sum","Up after Sum"))


mat <- countsHENorm[na.omit(int_up_down[1:7]),grepl("microglia|oligodendrocytes",cells$level1class)]
mat <- log(mat +1)
dfVar <- as.data.frame(colData(dds_sbst)[,"level1class"])
pheatmap(mat, annotation_col=dfVar,cluster_cols = FALSE)

mat_al <- countsHENorm_al[na.omit(int_up_down[1:7]),grepl("microglia|oligodendrocytes",cells$level1class)]
mat_al <- log(mat_al +1)
dfVar_al <- as.data.frame(colData(dds_sbst_al)[,"level1class"])
pheatmap(mat_al, annotation_col=dfVar_al,cluster_cols=FALSE)
