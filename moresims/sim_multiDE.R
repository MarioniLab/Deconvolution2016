# This checks the brittleness of the summation approach in the presence of DE genes.

source("functions.R")
dir.out <- "results_multiDE"
dir.create(dir.out, showWarnings=FALSE)
fout <- file.path(dir.out, "summary.txt")
is.there <- FALSE 

for (mode in c("UMI", "read")) { 
    for (overlap in c(TRUE, FALSE)) {
        cur.dir <- file.path(dir.out, paste0(mode, "-", overlap))
        dir.create(cur.dir, showWarnings=FALSE)
        ncells <- 1000

        for (genes.up in c(0, 1000, 2000)) {
            for (genes.down in c(0, 1000, 2000)) {
                if (genes.up==0 && genes.down==0) { next; }
                for (prop in c(0.1, 0.2, 0.3)) {
                    for (effect in c(2, 5, 10)) {

###### BROKEN INDENTING ######

    stub <- file.path(cur.dir, sprintf("u%i-d%i-p%.1f-e%i", genes.up, genes.down, prop, effect))
    param <- generateRawMeans(ncells=ncells, mode=mode)
    truth <- param$sf

    # Adding the requested DE (second population has reverse profile from the first population).
    mus <- param$fitted

    de.cell1 <- seq_len(ncells*prop)
    goes.up1 <- seq_len(genes.up)
    goes.down1 <- seq_len(genes.down) + genes.up
    mus[goes.up1,de.cell1] <- mus[goes.up1,de.cell1] * effect
    mus[goes.down1,de.cell1] <- 0

    de.cell2 <- seq_len(ncells*prop) + max(de.cell1)
    goes.up2 <- seq_len(genes.down) 
    if (!overlap) { # Do the DE genes in the second population overlap with the first?
        goes.up2 <- goes.up2 + max(goes.down1)
    }
    goes.down2 <- seq_len(genes.up) + max(goes.up2)
    mus[goes.up2,de.cell2] <- mus[goes.up2,de.cell2] * effect
    mus[goes.down2,de.cell2] <- 0

    de.cell <- c(de.cell1, de.cell2)
    counts <- sampleCounts(mus, param$dispersion)
    output <- runAllMethods(counts)

    # Generating an output plot for this simulation.
    collected <- vector("list", length(output))
    names(collected) <- names(output)
    col <- rep("black", ncells)
    col[de.cell1] <- "orange"
    col[de.cell2] <- "dodgerblue"

    pdf(paste0(stub, ".pdf"))
    par(mar=c(5.1,5.1,4.1,1.1))
    for (x in names(output)) {
        collected[[x]] <- makeSFPlot(output[[x]], truth, de.cell, main=x, col=col)["DE"]
    }
    dev.off()

    # Writing summary statistics to table.
    stats <- do.call(cbind, collected)
    write.table(file=fout, data.frame(Mode=mode, Overlap=overlap, Ncells=ncells, Up=genes.up, Down=genes.down, 
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

