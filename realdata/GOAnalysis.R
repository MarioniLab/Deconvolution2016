require("topGO")

# ---- Data ----

input.dir<- "./DEresults/edgeR"
lib <- read.table(file.path(input.dir,"ZeiselLibSize.tsv.gz"), header=TRUE, row.names=1) 
TMM <- read.table(file.path(input.dir,"ZeiselTMM.tsv.gz"), header=TRUE, row.names=1) 
sf <- read.table(file.path(input.dir,"ZeiselDESeq.tsv.gz"), header=TRUE, row.names=1) 
decon <- read.table(file.path(input.dir,"ZeiselDeconvolution.tsv.gz"), header=TRUE, row.names=1) 

sig.lib <- lib[lib$FDR < 0.05,]
sig.TMM <- TMM[TMM$FDR < 0.05,]
sig.sf <- sf[sf$FDR < 0.05,]
sig.decon <- decon[decon$FDR < 0.05,]

# ---- GOanalysis ----

## testing unique DE genes
uq.lib <- rownames(sig.lib)[!rownames(sig.lib) %in% rownames(sig.decon)]
uq.TMM <- rownames(sig.TMM)[!rownames(sig.TMM) %in% rownames(sig.decon)]
uq.sf <- rownames(sig.sf)[!rownames(sig.sf) %in% rownames(sig.decon)]
uq.decon <- rownames(sig.decon)[!rownames(sig.decon) %in% rownames(sig.lib)]

set.uq.lib <- as.integer(rownames(decon) %in% uq.lib)
set.uq.TMM <- as.integer(rownames(decon) %in% uq.TMM)
set.uq.sf <- as.integer(rownames(decon) %in% uq.sf)
set.uq.decon <- as.integer(rownames(decon) %in% uq.decon)
names(set.uq.lib) <- names(set.uq.TMM) <- names(set.uq.sf) <- names(set.uq.decon) <- rownames(decon)

#prepare Data for topGO

list.uq <- list("Lib"=set.uq.lib,"TMM"=set.uq.TMM,"SF"=set.uq.sf,"Decon"=set.uq.decon)
ontologies <- c("BP","CC","MF")
output.dir <- "GOresults"
dir.create(output.dir)

for (i in seq_along(list.uq)) {
    alg <- factor(list.uq[[i]])

    for (x in ontologies) {
        GO.data <- new("topGOdata", description="Lib GO",ontology=x, allGenes=alg, 
                      annot=annFUN.org, mapping="org.Mm.eg.db", nodeSize=10, ID="symbol")
        result.classic <- runTest(GO.data, algorithm="classic", statistic="Fisher")
        output <- GenTable(GO.data, Fisher.classic=result.classic, orderBy="Fisher.classic", topNodes=200, numChar=10000)
        write.table(output, file.path(output.dir, paste0("GOoutputUQin", names(list.uq)[i], "_", x, ".tsv")), quote=FALSE, row.names=FALSE, sep="\t")
    }
}
