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

scaled_factors <- data.frame("DESeq" = szf_HE / median(szf_HE),"TMM" = tmm/median(tmm), "Library size" = libfactor/median(libfactor), "Deconvolution" = szf_alClust/median(szf_alClust),check.names=FALSE)
#Comparison of SF Distribution

cairo_pdf("Zeisel_NormFactors.pdf")
par(mar=c(8.6,5.1,2.1,1.1))
boxplot(scaled_factors,log="y",ylab="Size factors",cex.axis=1.5,cex.lab=1.8,xaxt='n')
text(labels=colnames(scaled_factors),xpd=TRUE,cex=1.8, x=c(1.8,2.8,3.8,4.8)-0.7, y=0.1,srt=45, pos=2)
dev.off()
#SF vs Deconvoluted

cairo_pdf("Zeisel_SFvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,scaled_factors$Deconvolution,xlim=c(0.1,6),ylim=c(0.1,6),cex=0.6,pch=19,col="#00000073",xlab="DESeq",ylab="Deconvolution",log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")
dev.off()

#Libfactor vs Deconvoluted
cairo_pdf("Zeisel_LibvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors[,"Library size"],scaled_factors$Deconvolution,xlim=c(0.1,7),ylim=c(0.1,7),cex=0.6,pch=19,col="#00000073",xlab="Library size",ylab="Deconvolution",log="xy",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")
dev.off()

#TMM vs Deconvoluted
cairo_pdf("Zeisel_TMMvDeconv.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$TMM,scaled_factors$Deconvolution,xlim=c(0.1,7),ylim=c(0.1,7),cex=0.6,pch=19,col="#00000073",xlab="TMM",ylab="Deconvolution",cex.axis=1.5,cex.lab=1.8)
abline(0,1,col="dodgerblue")
dev.off()
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
cairo_pdf("Zeisel_ERCCvSF.pdf")
par(mar=c(5.1,5.1,4.1,1.1))
plot(scaled_factors$DESeq,erccfactor/median(erccfactor),ylab="Spike-in",xlab="DESeq",log="x",cex=0.6,pch=19,col="#00000073",cex.axis=1.5,cex.lab=1.8)
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


# ---- Diferential-Expression ----
## Subset Data for EdgeR
countssbst <- countsHEnoSpike[,grepl("microglia|oligodendrocytes",cells$level1class)]
cellssbst <- cells[grepl("microglia|oligodendrocytes",cells$level1class),]

szf_sbst<- szf_HE[grepl("microglia|oligodendrocytes",cells$level1class)]

szf_sbst_al<- szf_alClust[grepl("microglia|oligodendrocytes",cells$level1class)]

tmm_sbst<- tmm[grepl("microglia|oligodendrocytes",cells$level1class)]

lib_sbst <- libfactor[grepl("microglia|oligodendrocytes",cells$level1class)]

# Setup EdgeR
desgn <- model.matrix(~factor(cellssbst$level1class))
y <- DGEList(countssbst)

#Sizefactors
y.sf <- y
y.sf$samples$norm.factors <- szf_sbst/y.sf$samples$lib.size
y.sf <- estimateDisp(y.sf,desgn)
fit.sf <- glmFit(y.sf,desgn)
res.sf <- glmLRT(fit.sf)

#TMM
y.tmm <- y
y.tmm$samples$norm.factors <- tmm_sbst/y.tmm$samples$lib.size
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
y.al$samples$norm.factors <- szf_sbst_al/y.al$samples$lib.size
y.al <- estimateDisp(y.al,desgn)
fit.al <- glmFit(y.al,desgn)
res.al <- glmLRT(fit.al)

## Comparison

x.sf <- decideTestsDGE(res.sf)
x.lib <- decideTestsDGE(res.lib)
x.tmm <- decideTestsDGE(res.tmm)
x.al <- decideTestsDGE(res.al)


out.file <- "Zeisel_output.txt"
write.table(file=out.file, data.frame(Method="Size factor", Total=sum(x.sf!=0), Down=sum(x.sf<0), Up=sum(x.sf>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(file=out.file, data.frame(Method="TMM", Total=sum(x.tmm!=0), Down=sum(x.tmm<0), Up=sum(x.tmm>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="Library size", Total=sum(x.lib!=0), Down=sum(x.lib<0), Up=sum(x.lib>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="Deconvolution", Total=sum(x.al!=0), Down=sum(x.al<0), Up=sum(x.al>0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

write.table(file=out.file, data.frame(Method="vs SF", Total=sum(x.sf!=0 & x.al!=0), Down=sum(x.sf<0 & x.al < 0), Up=sum(x.sf>0 & x.al > 0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="vs TMM", Total=sum(x.tmm!=0 & x.al!=0), Down=sum(x.tmm<0 & x.al < 0), Up=sum(x.tmm>0 & x.al > 0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(file=out.file, data.frame(Method="vs lib", Total=sum(x.lib!=0 & x.al!=0), Down=sum(x.lib<0 & x.al < 0), Up=sum(x.lib>0 & x.al > 0)), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
