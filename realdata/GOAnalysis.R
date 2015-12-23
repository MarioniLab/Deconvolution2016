require("biomaRt")
require("topGO")

# ---- Data ----

input.dir<- "./DEresults"
lib <- read.table(file.path(input.dir,"ZeiselLibSize.tsv.gz"), header=TRUE, row.names=1) 
TMM <- read.table(file.path(input.dir,"ZeiselTMM.tsv.gz"), header=TRUE, row.names=1) 
sf <- read.table(file.path(input.dir,"ZeiselDESeq.tsv.gz"), header=TRUE, row.names=1) 
decon <- read.table(file.path(input.dir,"ZeiselDeconvolution.tsv.gz"), header=TRUE, row.names=1) 

sig.lib <- lib[lib$FDR < 0.05,]
sig.TMM <- TMM[TMM$FDR < 0.05,]
sig.sf <- sf[sf$FDR < 0.05,]
sig.decon <- decon[decon$FDR < 0.05,]


# ---- biomart ----

ensembl <- useMart(biomart= "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", host ="www.ensembl.org")
mm.genes <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),mart=ensembl)

geneConv <- mm.genes[mm.genes$external_gene_name %in% rownames(decon),]
geneConv <- mm.genes[match(rownames(decon),mm.genes$external_gene_name),]

# ---- GOanalysis ----

uq.lib <- rownames(sig.decon)[!rownames(sig.decon) %in% rownames(sig.lib)]
uq.TMM <- rownames(sig.decon)[!rownames(sig.decon) %in% rownames(sig.TMM)]
uq.sf <- rownames(sig.decon)[!rownames(sig.decon) %in% rownames(sig.sf)]

set.uq.lib <- as.integer(rownames(decon) %in% uq.lib)
set.uq.TMM <- as.integer(rownames(decon) %in% uq.TMM)
set.uq.sf <- as.integer(rownames(decon) %in% uq.sf)

list.uq <- list("Lib"=set.uq.lib,"TMM"=set.uq.TMM,"SF"=set.uq.sf)
ontologies <- c("BP","CC","MF")

output <- list("Lib"=as.list(ontologies),"TMM"=as.list(ontologies),"DESeq"=as.list(ontologies))
output.dir <- "GOresults"
dir.create(output.dir)


for (i in 1:length(list.uq)) {

    alg <- factor(list.uq[[i]])
    names(alg) <- geneConv$ensembl_gene_id
    
    for (x in ontologies) {

        GO.data <- new("topGOdata",
                      description="Lib GO",ontology=x,
                      allGenes=alg, 
                      annot=annFUN.org, mapping="org.Mm.eg.db",
                      nodeSize=5,ID="ensembl")

        result.classic <- runTest(GO.data, algorithm="classic", statistic="Fisher" )

        output[[i]][[x]] <- GenTable(GO.data, 
            Fisher.classic=result.classic,
            orderBy="Fisher.classic", topNodes = 100)
    }
    output.print <- data.frame(output[i])[,-1:-3]
    write.table(output.print,paste(output.dir,"/GOoutput",names(list.uq)[i],".txt",sep=""))
}
