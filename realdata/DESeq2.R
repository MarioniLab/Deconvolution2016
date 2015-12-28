require(DESeq2)
de.dir <- "DEresults/DESeq2"
dir.create(de.dir, recursive=TRUE)

for (x in c("Zeisel", "Klein")) {
    cur.data <- readRDS(sprintf("%sData.rds", x))
    counts <- cur.data$Counts
    size.facs <- cur.data$SF
    grouping <- cur.data$Cells

    if (x=="Zeisel") { 
        # Using only microglia
        of.interest <- grouping$level1class %in% c("pyramidal CA1", "oligodendrocytes")
        counts <- counts[,of.interest]
        size.facs <- size.facs[of.interest,]
        grouping <- grouping[of.interest,]
        y <- DESeqDataSetFromMatrix(countData = counts, colData = grouping, design = ~ level1class)
    } else {
        y <- DESeqDataSetFromMatrix(countData = counts, colData = grouping, design = ~ Time)
    }


    # Using the precomputed size factors.
    all.methods <- c("DESeq", "TMM", "LibSize", "Deconvolution")
    stopifnot(all(all.methods%in%colnames(size.facs)))
    for (method in all.methods) {

        colData(y)$sizeFactor <- size.facs[[method]]
        y <- estimateDispersions(y)
        y <- nbinomWaldTest(y)
        reslts <- results(y, alpha= 0.05)

        fhandle <- gzfile(file.path(de.dir, paste0(x, method, ".tsv.gz")), open="wb")
        write.table(file=fhandle, reslts, sep="\t", quote=FALSE, col.names=NA)
        close(fhandle)

        significant <- reslts$padj < 0.05
        significant[is.na(significant)] <- FALSE
        log2folds <- reslts$log2FoldChange
        log2folds[is.na(log2folds)] <- 0
        out <- log2folds * significant  #Create vector with log2foldchanges, 0 if not significant

        if (method=="DESeq") { 
            x.deseq <- out
        } else if (method=="TMM") {
            x.tmm <- out
        } else if (method=="LibSize") {
            x.lib <- out
        } else if (method=="Deconvolution") {
            x.d <- out
        }
        gc()
    }

    out.file <- file.path(de.dir, sprintf("%s_number.txt", x))
    save2file <- function(..., first=FALSE) { write.table(file=out.file, data.frame(...), sep="\t", quote=FALSE, row.names=FALSE, col.names=first, append=!first) }
    save2file(Method="DESeq", Total=sum(x.deseq!=0), Down=sum(x.deseq<0), Up=sum(x.deseq>0), first=TRUE)
    save2file(Method="TMM", Total=sum(x.tmm!=0), Down=sum(x.tmm<0), Up=sum(x.tmm>0))
    save2file(Method="Library size", Total=sum(x.lib!=0), Down=sum(x.lib<0), Up=sum(x.lib>0))
    save2file(Method="Deconvolution", Total=sum(x.d!=0), Down=sum(x.d<0), Up=sum(x.d>0))  
    save2file(Method="shared with DESeq", Total=sum(x.deseq!=0 & x.d!=0), Down=sum(x.deseq<0 & x.d < 0), Up=sum(x.deseq>0 & x.d > 0))
    save2file(Method="shared with TMM", Total=sum(x.tmm!=0 & x.d!=0), Down=sum(x.tmm<0 & x.d < 0), Up=sum(x.tmm>0 & x.d > 0))
    save2file(Method="shared with library size", Total=sum(x.lib!=0 & x.d!=0), Down=sum(x.lib<0 & x.d < 0), Up=sum(x.lib>0 & x.d > 0))
}
