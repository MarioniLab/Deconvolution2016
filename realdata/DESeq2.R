require(DESeq2)

for (x in c("Zeisel","Klein")) { 
    counts <- read.csv(sprintf("%sCounts.csv", x), row.names=1)
    size.facs <- read.csv(sprintf("%sSF.csv", x))
    grouping <- read.csv(sprintf("%sCells.csv", x))

    if (x=="Zeisel") { 
        # Using only microglia
        of.interest <- grouping$level1class %in% c("microglia", "oligodendrocytes")
        counts <- counts[,of.interest]
        size.facs <- size.facs[of.interest,]
        grouping <- grouping[of.interest,]
        y <- DESeqDataSetFromMatrix(countData = counts, colData = grouping, design = ~ level1class)
    } else {
        y <- DESeqDataSetFromMatrix(countData = counts, colData = grouping, design = ~ Cells)
    }


    # Using the precomputed size factors.
    for (method in c("SF", "TMM", "lib", "Decon")) {

        colData(y)$sizeFactor <- size.facs[[method]]
        y <- estimateDispersions(y)
        y <- nbinomWaldTest(y)
        reslts <- results(y, alpha= 0.05)

        significant <- reslts$padj < 0.05
        significant[is.na(significant)] <- FALSE
        log2folds <- reslts$log2FoldChange
        log2folds[is.na(log2folds)] <- 0
        out <- log2folds * significant  #Create vector with log2foldchanges, 0 if not significant

        if (method=="SF") { 
            x.sf <- out
        } else if (method=="TMM") {
            x.tmm <- out
        } else if (method=="lib") {
            x.lib <- out
        } else {
            x.d <- out
        }
        gc()
    }

    out.file <- sprintf("%s_output.txt", x)
    save2file <- function(..., first=FALSE) { write.table(file=out.file, data.frame(...), sep="\t", quote=FALSE, row.names=FALSE, col.names=first, append=!first) }
    save2file(Method="DESeq", Total=sum(x.sf!=0), Down=sum(x.sf<0), Up=sum(x.sf>0), first=TRUE)
    save2file(Method="TMM", Total=sum(x.tmm!=0), Down=sum(x.tmm<0), Up=sum(x.tmm>0))
    save2file(Method="Library size", Total=sum(x.lib!=0), Down=sum(x.lib<0), Up=sum(x.lib>0))
    save2file(Method="Deconvolution", Total=sum(x.d!=0), Down=sum(x.d<0), Up=sum(x.d>0))  
    save2file(Method="shared with DESeq", Total=sum(x.sf!=0 & x.d!=0), Down=sum(x.sf<0 & x.d < 0), Up=sum(x.sf>0 & x.d > 0))
    save2file(Method="shared with TMM", Total=sum(x.tmm!=0 & x.d!=0), Down=sum(x.tmm<0 & x.d < 0), Up=sum(x.tmm>0 & x.d > 0))
    save2file(Method="shared with library size", Total=sum(x.lib!=0 & x.d!=0), Down=sum(x.lib<0 & x.d < 0), Up=sum(x.lib>0 & x.d > 0))
}
