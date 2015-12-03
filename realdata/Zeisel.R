# AK47 Normalization

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
szf_HE <- estimateSizeFactorsForMatrix(countsHEnoSpike,geoMeans = geoMeans_HE)

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

scaled_factors <- data.frame("Sizefactor" = szf_HE / median(szf_HE),"TMM" = tmm/median(tmm), "LibrarySize" = libfactor/median(libfactor), "Deconvoluted" = szf_alClust/median(szf_alClust))
#Comparison of SF Distribution
boxplot(scaled_factors,log="y",ylab="Scaled Normalizationfactor")

#SF vs Deconvoluted
plot(scaled_factors$Sizefactor,scaled_factors$Deconvoluted,xlim=c(0.1,6),ylim=c(0.1,6),cex=0.6,pch=19,col="#00000073",xlab="Sizefactor",ylab="Deconvoluted",log="xy")
abline(0,1,col="dodgerblue")

#Libfactor vs Deconvoluted
plot(scaled_factors$LibrarySize,scaled_factors$Deconvoluted,xlim=c(0.1,7),ylim=c(0.1,7),cex=0.6,pch=19,col="#00000073",xlab="Librarysize",ylab="Deconvoluted",log="xy")
abline(0,1,col="dodgerblue")

#TMM vs Deconvoluted
plot(scaled_factors$TMM,scaled_factors$Deconvoluted,xlim=c(0.1,7),ylim=c(0.1,7),cex=0.6,pch=19,col="#00000073",xlab="TMM",ylab="Deconvoluted")
abline(0,1,col="dodgerblue")

# ---- Normalize ----

# countsHENorm <- as.data.frame(t(t(countsHEnoSpike) / szf_HE))
# countsHENorm_alClust <- as.data.frame(t(t(countsHEnoSpike) / szf_alClust))

# ---- LibrarySize ----
# Librarysizes <- c(colSums(countsHE),colSums(countsHENorm),colSums(countsHENorm_alClust))
# Normalization <- c(rep("Raw",ncol(countsHE)),rep("SF",ncol(countsHENorm)),rep("SF_al",ncol(countsHENorm_alClust)))
# ddf <- data.frame(Librarysizes,Normalization)
# plt <- ggplot(ddf, aes(x=Librarysizes, fill=Normalization, group=Normalization))
# plot(plt + geom_density(alpha=0.3) + scale_x_log10())
# ---- Spike-Ins ----

# plot(scaled_factors$Sizefactor,erccfactor/median(erccfactor),ylab="Spike-In Sizefactor",xlab="Sizefactors",log="x",cex=0.6,pch=19,col="#00000073")
# plot(scaled_factors$TMM,erccfactor/median(erccfactor),ylab="Spike-In Sizefactor",xlab="TMM",log="x",cex=0.6,pch=19,col="#00000073")
# plot(scaled_factors$Deconvoluted,erccfactor/median(erccfactor),ylab="Spike-In Sizefactor",xlab="Deconvoluted",log="x",cex=0.6,pch=19,col="#00000073")
# plot(scaled_factors$LibrarySize,erccfactor/median(erccfactor),ylab="Spike-In Sizefactor",xlab="LibrarySize",log="x",cex=0.6,pch=19,col="#00000073")

## Compute conversionFactor
# betafactor <- apply(countsHE[featuresHE$spike,],2, function(x) lm(x ~ molecules[molecules$ERCC_ID %in% rownames(featuresHE)[featuresHE$spike],"molecules_in_each_chamber"])$coefficients[2])

## Sizefactors dont correspond to conversion factors
# plot(szf_HE,betafactor,cex=0.6,pch=19,col="#00000073",xlab="Sizefactors",ylab="Conversion Factors",main="Sizefactors are not assoicated with Conversion Factors")

# ---- Effect-on-HVG ---

##Order of Variability
# plot(order(dm),order(dmNorm),cex=0.5,pch=19)
# plot(order(dm),order(dmNorm_al),cex=0.5,pch=19)
# plot(order(dmNorm),order(dmNorm_al),cex=0.5,pch=19)
# 
# featuresHE_Ord<- featuresHE[order(-dm),]
# featuresNorm_Ord<- featuresNorm[order(-dmNorm),]
# featuresNorm_alOrd <- featuresNorm_al[order(-dmNorm_al),]
# 
# featuresHE$hvg <- rownames(featuresHE) %in% rownames(featuresHE_Ord)[1:1000]
# featuresNorm$hvg <- rownames(featuresNorm) %in% rownames(featuresNorm_Ord)[1:1000]
# featuresNorm_al$hvg <- rownames(featuresNorm_al) %in% rownames(featuresNorm_alOrd)[1:1000]
# 
# venHE <- which(featuresHE$hvg)
# venNorm <- which(featuresNorm$hvg)
# venNorm_al <- which(featuresNorm_al$hvg)
# 
# venn.diagram(list("Non-Norm" = venHE,"Sizefacotrs" = venNorm, "Summation" = venNorm_al), fill = c("red","green","blue"),
#              alpha= c(0.5,0.5,0.5), cex = 1 , cat.fontface = 4, fontfamily = 3,imagetype="png",
#              filename = "VenDiagram2.png")

# ---- Diferential-Expression ----

sbst <- countsHEnoSpike[,grepl("microglia|oligodendrocytes",cells$level1class)]
cellssbst <- cells[grepl("microglia|oligodendrocytes",cells$level1class),]
cellssbst$sizeFactor <- szf_HE[grepl("microglia|oligodendrocytes",cells$level1class)]

cellssbst_al <- cells[grepl("microglia|oligodendrocytes",cells$level1class),]
cellssbst_al$sizeFactor <- szf_alClust[grepl("microglia|oligodendrocytes",cells$level1class)]

cellssbst_tmm <- cells[grepl("microglia|oligodendrocytes",cells$level1class),]
cellssbst_tmm$sizeFactor <- tmm[grepl("microglia|oligodendrocytes",cells$level1class)]

cellssbst_lib<- cells[grepl("microglia|oligodendrocytes",cells$level1class),]
cellssbst_lib$sizeFactor <- libfactor[grepl("microglia|oligodendrocytes",cells$level1class)]

#DESeq for SF
dds_sbst <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst, design =~ level1class)
dds_sbst <- estimateDispersions(dds_sbst)
dds_sbst <- nbinomWaldTest(dds_sbst)
rslt_sbst <- results(dds_sbst)
write.csv(as.data.frame(rslt_sbst),"Zeisel_Result_SF.csv")

#DESeq for Deconvoluted
dds_sbst_al <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst_al, design =~ level1class)
dds_sbst_al <- estimateDispersions(dds_sbst_al)
dds_sbst_al <- nbinomWaldTest(dds_sbst_al)
rslt_sbst_al <- results(dds_sbst_al)
write.csv(as.data.frame(rslt_sbst_al),"Zeisel_Result_Deconv.csv")

#DESeq for TMM
dds_sbst_tmm <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst_tmm, design =~ level1class)
dds_sbst_tmm <- estimateDispersions(dds_sbst_tmm)
dds_sbst_tmm <- nbinomWaldTest(dds_sbst_tmm)
rslt_sbst_tmm <- results(dds_sbst_tmm)
write.csv(as.data.frame(rslt_sbst_tmm),"Zeisel_Result_TMM.csv")

#DESeq for Libfactor
dds_sbst_lib <- DESeqDataSetFromMatrix(countData = sbst, colData =cellssbst_lib, design =~ level1class)
dds_sbst_lib <- estimateDispersions(dds_sbst_lib)
dds_sbst_lib <- nbinomWaldTest(dds_sbst_lib)
rslt_sbst_lib <- results(dds_sbst_lib)
write.csv(as.data.frame(rslt_sbst_lib),"Zeisel_Result_Lib.csv")

# up_sbst <- rownames(rslt_sbst)[rslt_sbst$log2FoldChange > 0 & rslt_sbst$padj < 0.1]
# up_sbst_al <- rownames(rslt_sbst_al)[rslt_sbst_al$log2FoldChange > 0 & rslt_sbst_al$padj < 0.1]
# 
# down_sbst <- rownames(rslt_sbst)[rslt_sbst$log2FoldChange < 0 & rslt_sbst$padj < 0.1]
# down_sbst_al <- rownames(rslt_sbst_al)[rslt_sbst_al$log2FoldChange < 0 & rslt_sbst$padj < 0.1]
# 
# int_up <- intersect(up_sbst,up_sbst_al)
# int_down <- intersect(down_sbst, down_sbst_al)
# 
# int_up_down <- intersect(up_sbst,down_sbst_al)
# int_down_up <- intersect(down_sbst,up_sbst_al)
# 
# barplot(c(length(up_sbst),length(up_sbst_al),length(int_up),length(down_sbst),length(down_sbst_al),length(int_down),length(int_up_down),length(int_down_up)),width=0.5,names.arg=c("SF","SF_AS","Overlap","SF","SF_AS","Overlap","Overlap","Overlap"),col=c("#FB6A4A","#FB6A4A","#EF3B2C","#6BAED6","#6BAED6","#2171B5","seagreen2","tomato2"),cex.names=0.7,space=c(0.1,0.1,0.1,0.5,0.1,0.1,0.5,0.1),main="Oligodendrocyte v. Microglia")
# legend("topright",fill=c("#FB6A4A","#6BAED6","seagreen2","tomato2"),legend=c("Upregulated","Downregulated","Down after Sum","Up after Sum"))

