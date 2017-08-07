# This checks the brittleness of the summation approach in the presence of DE genes.

source("functions.R")
dir.out <- "results_noDE"
dir.create(dir.out, showWarnings=FALSE)
fout <- file.path(dir.out, "summary.txt")
is.there <- FALSE 

for (mode in c("UMI", "read")) { 
    dir.create(file.path(dir.out, mode), showWarnings=FALSE)
    for (ncells in c(100, 200, 500, 1000)) { 
        stub <- file.path(dir.out, mode, ncells)

###### BROKEN INDENTING ######

    param <- generateRawMeans(ncells=ncells, mode=mode)
    truth <- param$sf
    counts <- sampleCounts(param$fitted, param$dispersion)
    output <- runAllMethods(counts)

    # Generating an output plot for this simulation.
    pdf(paste0(stub, ".pdf"))
    par(mar=c(5.1,5.1,4.1,1.1))
    collected <- vector("list", length(output))
    names(collected) <- names(output)
    for (x in names(output)) {
        out <- makeSFPlot(output[[x]], truth, main=x)
        collected[[x]] <- out["non-DE"]
    }
    dev.off()

    # Writing summary statistics to table.
    stats <- do.call(cbind, collected)
    write.table(file=fout, data.frame(Mode=mode, Ncells=ncells, format(stats, digits=3)),
                append=is.there, quote=FALSE, sep="\t", col.names=!is.there, row.names=FALSE)
    is.there <- TRUE

###### END INDENTING ######

    }
}

