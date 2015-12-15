require(edgeR)
dir.create("DEresults")

for (x in c("Zeisel", "Klein")) { 
    counts <- read.csv(sprintf("%sCounts.csv", x), row.names=1)
    size.facs <- read.csv(sprintf("%sSF.csv", x))
    grouping <- read.csv(sprintf("%sCells.csv", x))

    if (x=="Zeisel") { 
        # Using only microglia
        of.interest <- grouping$level1class %in% c("pyramidal CA1", "oligodendrocytes")
        counts <- counts[,of.interest]
        size.facs <- size.facs[of.interest,]
        grouping <- grouping[of.interest,]
        design <- model.matrix(~factor(grouping$level1class))
    } else {
        design <- model.matrix(~factor(grouping$Cells))        
    }

    # Using the precomputed size factors.
    for (method in c("SF", "TMM", "lib", "Decon")) {
        y <- DGEList(counts)
        y$samples$norm.factors <- size.facs[[method]]/y$samples$lib.size
        y <- estimateDisp(y, design, prior.df=0, trend='none')
        fit <- glmFit(y, design)
        res <- glmTreat(fit, lfc=1)

        xxx <- topTags(res, n=Inf, sort.by="none")
        fhandle <- gzfile(file.path("DEresults", paste0(x, method, ".tsv.gz")), open="wb")
        write.table(file=fhandle, xxx$table, sep="\t", quote=FALSE, col.names=NA)
        close(fhandle)

        out <- decideTestsDGE(res)
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
