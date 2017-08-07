##################################################
# Setting up a comparator function.

library(viridis)
inspector <- function(blah, split.levels, to.compare) {
    all.f <- vector("list", length(split.levels))
    for (fi in seq_along(split.levels)) {
        all.f[[fi]] <- paste(split.levels[fi], "=", blah[,split.levels[fi]])
    }
    f <- do.call(paste, c(all.f, sep=", "))
    by.f <- split(blah, f)

    for (x in names(by.f)) {
        current <- by.f[[x]]
        left <- current[[to.compare[1]]]*100
        right <- current[[to.compare[2]]]*100
        max.pt <- max(left, right)

        # Determine balance by colour, number of DE genes by size.
        if ("Up" %in% names(current)) { 
            balance <- factor(current$Up - current$Down)
            col <- viridis(nlevels(balance))[balance]
            num.de <- factor(current$Up + current$Down)
            effect.size <- factor(current$Effect)
        } else {
            col <- "black"
            num.de <- 1
            effect.size <- 1
        }

        plot(left, right, xlim=c(0, max.pt), ylim=c(0, max.pt), 
             col=col, pch=as.integer(effect.size), cex=as.integer(num.de),
             xlab=paste(to.compare[1], "error (%)"), 
             ylab=paste(to.compare[2], "error (%)"), main=x)
        abline(0, 1, col="red", lty=2)
    }
    return(invisible(NULL))
}

##################################################
# This visually compares the error of library size normalization to deconvolution (w/o clustering).

blah <- read.table("results_biDE/summary.txt", header=TRUE, stringsAsFactors=FALSE)
pdf("results_biDE/summary_deconv_vs_lib.pdf")
inspector(blah, c("Mode", "Ncells"), c("Lib", "Deconv"))
dev.off()

pdf("results_biDE/summary_deconv_vs_clust.pdf")
inspector(blah, c("Mode", "Ncells"), c("Deconv", "Deconv.clust"))
dev.off()

blah <- read.table("results_multiDE/summary.txt", header=TRUE, stringsAsFactors=FALSE)
pdf("results_multiDE/summary_deconv_vs_lib.pdf")
inspector(blah, c("Mode", "Overlap"), c("Lib", "Deconv"))
dev.off()

pdf("results_multiDE/summary_deconv_vs_clust.pdf")
inspector(blah, c("Mode", "Overlap"), c("Deconv", "Deconv.clust"))
dev.off()
