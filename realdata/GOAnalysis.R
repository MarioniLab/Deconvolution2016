require("biomaRt")
require("topGO")

# ---- Data ----

input.dir<- "./DEresults"
lib <- read.table(file.path(input.dir,"Zeisellib.tsv.gz"), header=TRUE, row.names=1) 
TMM <- read.table(file.path(input.dir,"ZeiselTMM.tsv.gz"), header=TRUE, row.names=1) 
sf <- read.table(file.path(input.dir,"ZeiselSF.tsv.gz"), header=TRUE, row.names=1) 
decon <- read.table(file.path(input.dir,"ZeiselDecon.tsv.gz"), header=TRUE, row.names=1) 

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

## testing unique DE genes
uq.lib <- rownames(sig.lib)[!rownames(sig.lib) %in% rownames(sig.decon)]
uq.TMM <- rownames(sig.TMM)[!rownames(sig.TMM) %in% rownames(sig.decon)]
uq.sf <- rownames(sig.sf)[!rownames(sig.sf) %in% rownames(sig.decon)]
uq.decon <- rownames(sig.decon)[!rownames(sig.decon) %in% rownames(sig.lib)]

set.uq.lib <- as.integer(rownames(decon) %in% uq.lib)
set.uq.TMM <- as.integer(rownames(decon) %in% uq.TMM)
set.uq.sf <- as.integer(rownames(decon) %in% uq.sf)
set.uq.decon <- as.integer(rownames(decon) %in% uq.decon)
# 
#testint OG up (that is down) regulated genes
# up.sf <- rownames(sig.sf)[sig.sf$logFC < 0]
# up.TMM <- rownames(sig.TMM)[sig.TMM$logFC < 0]
# up.lib <- rownames(sig.lib)[sig.lib$logFC < 0]
# up.decon <- rownames(sig.decon)[sig.decon$logFC < 0]
# 
##Oligodendrocyte marker genes
# marker <- c("Olig1","Olig2","Olig3","Cnp","Cntnap2","Mag","Mog","Ascl1","Rtn4r","Omg","Mbp","Pdgfra")
# 
# set.up.sf <- as.integer(rownames(decon) %in% up.sf)
# set.up.TMM <- as.integer(rownames(decon) %in% up.TMM)
# set.up.lib <- as.integer(rownames(decon) %in% up.lib)
# set.up.decon <- as.integer(rownames(decon) %in% up.decon)
# 
#prepare Data for topGO

list.uq <- list("Lib"=set.uq.lib,"TMM"=set.uq.TMM,"SF"=set.uq.sf,"Decon"=set.uq.decon)
ontologies <- c("BP","CC","MF")

output <- list("Lib"=as.list(ontologies),"TMM"=as.list(ontologies),"DESeq"=as.list(ontologies),"Decon"=as.list(ontologies))
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
                      nodeSize=10,ID="ensembl")

        result.classic <- runTest(GO.data, algorithm="classic", statistic="Fisher" )

        output[[i]][[x]] <- GenTable(GO.data, 
            Fisher.classic=result.classic,
            orderBy="Fisher.classic", topNodes = 200)
    }
    output.print <- data.frame(output[i])[,-1:-3]
    write.csv(output.print[,2:4],paste(output.dir,"/GOoutputUQin",names(list.uq)[i],".csv",sep=""))
}
