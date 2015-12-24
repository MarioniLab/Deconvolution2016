require(edgeR)
de.dir <- "DEresults"
dir.create(de.dir)

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
        design <- model.matrix(~factor(grouping$level1class))
    } else {
        design <- model.matrix(~factor(grouping$Time))  
    }

    # Using the precomputed size factors.
    all.methods <- c("DESeq", "TMM", "LibSize", "Deconvolution")
    stopifnot(all(all.methods%in%colnames(size.facs)))
    for (method in all.methods) {
        y <- DGEList(counts)
        y$samples$norm.factors <- size.facs[[method]]/y$samples$lib.size
        y <- estimateDisp(y, design, prior.df=0, trend='none')
        fit <- glmFit(y, design)
        res <- glmTreat(fit, lfc=1)

        xxx <- topTags(res, n=Inf, sort.by="none")
        fhandle <- gzfile(file.path(de.dir, paste0(x, method, ".tsv.gz")), open="wb")
        write.table(file=fhandle, xxx$table, sep="\t", quote=FALSE, col.names=NA)
        close(fhandle)

        out <- decideTestsDGE(res)
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

    out.file <- sprintf("%s_edgeR_number.txt", x)
    save2file <- function(..., first=FALSE) { write.table(file=out.file, data.frame(...), sep="\t", quote=FALSE, row.names=FALSE, col.names=first, append=!first) }
    save2file(Method="DESeq", Total=sum(x.deseq!=0), Down=sum(x.deseq<0), Up=sum(x.deseq>0), first=TRUE)
    save2file(Method="TMM", Total=sum(x.tmm!=0), Down=sum(x.tmm<0), Up=sum(x.tmm>0))
    save2file(Method="Library size", Total=sum(x.lib!=0), Down=sum(x.lib<0), Up=sum(x.lib>0))
    save2file(Method="Deconvolution", Total=sum(x.d!=0), Down=sum(x.d<0), Up=sum(x.d>0))  
    save2file(Method="shared with DESeq", Total=sum(x.deseq!=0 & x.d!=0), Down=sum(x.deseq<0 & x.d < 0), Up=sum(x.deseq>0 & x.d > 0))
    save2file(Method="shared with TMM", Total=sum(x.tmm!=0 & x.d!=0), Down=sum(x.tmm<0 & x.d < 0), Up=sum(x.tmm>0 & x.d > 0))
    save2file(Method="shared with library size", Total=sum(x.lib!=0 & x.d!=0), Down=sum(x.lib<0 & x.d < 0), Up=sum(x.lib>0 & x.d > 0))

    # Computing the rankings.
    
    out.deseq <- read.table(file.path(de.dir, paste0(x, all.methods[1], ".tsv.gz")), header=TRUE, row.names=1) 
    out.tmm <- read.table(file.path(de.dir, paste0(x, all.methods[2], ".tsv.gz")), header=TRUE, row.names=1) 
    out.lib <- read.table(file.path(de.dir, paste0(x, all.methods[3], ".tsv.gz")), header=TRUE, row.names=1) 
    out.decon <- read.table(file.path(de.dir, paste0(x, all.methods[4], ".tsv.gz")), header=TRUE, row.names=1) 
    
    r.deseq <- rank(out.deseq$PValue, ties="first") # WARNING: p-value not accurate for < 100 genes, as they're all zero!
    r.tmm <- rank(out.tmm$PValue, ties="first")
    r.lib <- rank(out.lib$PValue, ties="first")
    r.decon <- rank(out.decon$PValue, ties="first")

    out.file <- sprintf("%s_edgeR_ranking.txt", x)
    isfirst <- TRUE
    comp <- function(a, b, top) { sprintf("%.2f", sum(a<=top & b<=top)/top) }
    for (top in c(100, 500, 2000)) {
        write.table(file=out.file, data.frame(Top=top, SF=comp(r.deseq, r.decon, top), TMM=comp(r.tmm, r.decon, top), lib=comp(r.lib, r.decon, top)),
                    append=!isfirst, col.names=isfirst, row.names=FALSE, quote=FALSE, sep="\t")
        isfirst <- FALSE
    }
}
