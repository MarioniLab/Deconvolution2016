# This checks the brittleness of the summation approach in the presence of DE genes.

source("functions.R")
dir.out <- "results_noDE"
dir.create(dir.out, showWarnings=FALSE)

for (mode in c("UMI", "read")) { 
    dir.create(file.path(dir.out, mode), showWarnings=FALSE)
    for (ncells in c(100, 200, 500, 1000)) { 
        stub <- file.path(dir.out, mode, ncells)

###### BROKEN INDENTING ######

    param <- generateRawMeans(ncells=ncells, mode=mode)
    truth <- param$sf
    counts <- sampleCounts(param$mu, param$phi)
    output <- runAllMethods(counts)

    # Generating an output plot for this simulation.
    pdf(paste0(stub, ".pdf"))
    par(mar=c(5.1,5.1,4.1,1.1))
    returned <- numeric(length(output))
    names(returned) <- names(output)
    for (x in names(output)) {
        r2 <- makeSFPlot(output[[x]], truth, main=x)
        returned[[x]] <- r2
    }
    dev.off()

    # Printing out some stats.
    write.table(data.frame(Method=names(returned), Ncells=ncells, R2=unlist(returned)),
                file=paste0(stub, ".tab"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)                               

###### END INDENTING ######

    }
}

