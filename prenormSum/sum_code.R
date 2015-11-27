# This function sorts cells by their library sizes, and generates an ordering vector.

.generateSphere <- function(lib.sizes) {
    nlibs <- length(lib.sizes)
    o <- order(lib.sizes)
    even <- seq(2,nlibs,2)
    odd <- seq(1,nlibs,2)
    out <- c(o[odd], rev(o[even]))
    c(out, out)
}

# This contains the function that performs normalization on the summed counts.

normalizeBySums <- function(counts, sizes=c(20, 40, 60, 80, 100), positive=FALSE) {
    ncells <- ncol(counts)
    if (any(sizes >= ncells)) { 
        stop("not enough cells for specified 'sizes'") 
    }

    # Computing the necessary statistics.
    lib.sizes <- colSums(counts)
    sphere <- .generateSphere(lib.sizes)
    exprs <- t(t(counts)/lib.sizes)
    ave.cell <- rowMeans(exprs)

    # Using our summation approach.
    design <- list()
    output <- list()

    for (size in sizes) {
        all.mat <- matrix(0L, ncells, ncells)
        out.nf <- numeric(ncells)

        for (it in seq_len(ncells)){
            chosen <- sphere[it:(it+size-1L)] 
            all.mat[it,chosen] <- 1L
            cur.combined <- rowSums(exprs[,chosen])
            out.nf[it] <- median(cur.combined/ave.cell)
        }

        design <- c(design, list(all.mat))
        output <- c(output, list(out.nf))
    }

    design <- do.call(rbind, design)
    output <- unlist(output)

    # Adding extra equations to guarantee solvability (downweighted).
    weights <- rep(c(1, 0.00001), c(nrow(design), ncells))
    design <- rbind(design, diag(ncells))
    output <- c(output, apply(exprs/ave.cell, 2, median))
    root.weights <- sqrt(weights)

    if (positive) { 
        fitted <- limSolve::lsei(A=design*root.weights, B=output*root.weights, G=diag(ncells), H=numeric(ncells), type=2)
        final.nf <- fitted$X
    } else {
        final.nf <- solve(qr(design * root.weights), output * root.weights)
    }

    # Returning size factors, rather than normalization factors.
    final.sf <- final.nf * lib.sizes
    final.sf <- final.sf/mean(final.sf[final.sf>0])
    return(final.sf)
}



