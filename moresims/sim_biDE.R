# This checks the brittleness of the summation approach in the presence of DE genes.

source("functions.R")
dir.out <- "results_biDE"
dir.create(dir.out, showWarnings=FALSE)
fout <- file.path(dir.out, "summary.txt")
is.there <- FALSE 

for (mode in c("UMI", "read")) { 
    for (ncells in c(100, 200, 500, 1000)) { 
        cur.dir <- file.path(dir.out, paste0(mode, "-", ncells))
        dir.create(cur.dir, showWarnings=FALSE)

        for (genes.up in c(0, 1000, 2000)) {
            for (genes.down in c(0, 1000, 2000)) {
                if (genes.up==0 && genes.down==0) { next; }
                for (prop in c(0.1, 0.2, 0.5)) {
                    for (effect in c(2, 5, 10)) {

###### BROKEN INDENTING ######

    stub <- file.path(cur.dir, sprintf("u%i-d%i-p%.1f-e%i", genes.up, genes.down, prop, effect))
    param <- generateRawMeans(ncells=ncells, mode=mode)
    truth <- param$sf

    # Adding the requested DE.
    mus <- param$fitted
    de.cell <- seq_len(ncells*prop)
    goes.up <- seq_len(genes.up)
    goes.down <- seq_len(genes.down) + genes.up

    mus[goes.up,de.cell] <- mus[goes.up,de.cell] * effect
    mus[goes.down,de.cell] <- 0

    counts <- sampleCounts(mus, param$dispersion)
    output <- runAllMethods(counts)

    # Generating an output plot for this simulation.
    collected <- vector("list", length(output))
    names(collected) <- names(output)
    col <- rep("black", ncells)
    col[de.cell] <- "red"

    pdf(paste0(stub, ".pdf"))
    par(mar=c(5.1,5.1,4.1,1.1))
    for (x in names(output)) {
        collected[[x]] <- makeSFPlot(output[[x]], truth, de.cell, main=x, col=col)["DE"]
    }
    dev.off()

    # Writing summary statistics to table.
    stats <- do.call(cbind, collected)
    write.table(file=fout, data.frame(Mode=mode, Ncells=ncells, Up=genes.up, Down=genes.down, 
                                      Prop=prop, Effect=effect, format(stats, digits=3)),
                append=is.there, quote=FALSE, sep="\t", col.names=!is.there, row.names=FALSE)
    is.there <- TRUE

###### END INDENTING ######

                    }
                }
            }
        }
    }       
}

