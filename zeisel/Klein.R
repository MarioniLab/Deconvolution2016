### Normalization without Spike-Ins
source("functions.R")

esd0 <- read.csv("../data/Klein/GSM1599494_ES_d0_main.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  
esd2 <- read.csv("../data/Klein/GSM1599497_ES_d2_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  
esd4 <- read.csv("../data/Klein/GSM1599498_ES_d4_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  
esd7 <- read.csv("../data/Klein/GSM1599499_ES_d7_LIFminus.csv.bz2", header = FALSE, stringsAsFactors= FALSE)  

colnames(esd0) <- c("feature",c(1:933))
colnames(esd2) <- c("feature",c(934:1236))
colnames(esd4) <- c("feature",c(1237:1919))
colnames(esd7) <- c("feature",c(1920:2717))

counts <- merge (esd0,esd2,by="feature",all=T)
counts <- merge (counts, esd4, by="feature", all=T)
counts <- merge (counts, esd7, by="feature", all=T)
rownames(counts) <- counts$feature
counts <- counts[,!(colnames(counts) %in% "feature")]

##Construct feature and cell data
celltype<- c(rep("d0",933),rep("d2", 303),rep("d4", 683),rep("d7" , 798))
cells <- cellCalc(counts = counts, celltype = celltype,exprmin =1) 
