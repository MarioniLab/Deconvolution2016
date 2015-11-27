# A stub for a normalization function for CPMs.

.normfun <- function(x1, x2) {
    return(median(x1/x2))
}

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

normalizeBySums <- function(counts, sizes=c(20, 40, 60, 80, 100), 
                            clusters=NULL, positive=FALSE, rescale=TRUE) {
    ncells <- ncol(counts)
    if (!is.null(clusters)) {
        if (ncells!=length(clusters)) { 
            stop("'counts' ncols is not equal to 'clusters' length")
        }
        indices <- split(seq_len(ncells), clusters)
    } else {
        indices <- list(seq_len(ncells))
    }

    sizes <- as.integer(sizes)
    if (anyDuplicated(sizes)) { 
        stop("'sizes' is not unique") 
    } 

    # Computing the necessary statistics.
    lib.sizes <- colSums(counts)
    exprs <- t(t(counts)/lib.sizes)
    clust.nf <- clust.profile <- clust.libsize <- list()

    # Computing normalization factors within each cluster first.    
    for (clust in seq_along(indices)) { 
        curdex <- indices[[clust]]
        cur.exprs <- exprs[,curdex,drop=FALSE]
        cur.libs <- lib.sizes[curdex]
        cur.cells <- length(curdex)

        # Checking cluster sizes
        if (any(sizes > cur.cells)) { 
            stop("not enough cells for specified 'sizes'") 
        } else if (any(sizes*2L > cur.cells)) {
            warning("number of cells should be at least twice that of the largest 'sizes'")
        }

        # Using our summation approach.
        sphere <- .generateSphere(cur.libs)
        ave.cell <- rowMeans(cur.exprs)
        design <- list()
        output <- list()

        for (size in sizes) {
            all.mat <- matrix(0L, cur.cells, cur.cells)
            out.nf <- numeric(cur.cells)

            for (it in seq_len(cur.cells)){
                chosen <- sphere[it:(it+size-1L)] 
                all.mat[it,chosen] <- 1L
                cur.combined <- rowSums(cur.exprs[,chosen])
                out.nf[it] <- .normfun(cur.combined, ave.cell)
            }

            design <- c(design, list(all.mat))
            output <- c(output, list(out.nf))
        }

        design <- do.call(rbind, design)
        output <- unlist(output)

        # Adding extra equations to guarantee solvability (downweighted).
        weights <- rep(c(1, 0.00001), c(nrow(design), cur.cells))
        design <- rbind(design, diag(cur.cells))
        output <- c(output, apply(cur.exprs/ave.cell, 2, median))
        root.weights <- sqrt(weights)

        if (positive) { 
            fitted <- limSolve::lsei(A=design*root.weights, B=output*root.weights, G=diag(cur.cells), H=numeric(cur.cells), type=2)
            final.nf <- fitted$X
        } else {
            final.nf <- solve(qr(design * root.weights), output * root.weights)
            if (any(final.nf < 0)) { stop("negative factor estimates, re-run with 'positive=TRUE'") }
        }

        # Adding per-cluster information.
        clust.nf[[clust]] <- final.nf
        clust.profile[[clust]] <- ave.cell
        clust.libsize[[clust]] <- mean(cur.libs)
    }

    # Adjusting size factors between clusters (using the cluster with the median per-cell library size as the reference).
    clust.libsize <- unlist(clust.libsize)
    ref.col <- which(rank(clust.libsize, ties.method="first")==as.integer(length(clust.libsize)/2)+1L)
    for (clust in seq_along(indices)) { 
        clust.nf[[clust]] <- clust.nf[[clust]] * .normfun(clust.profile[[clust]], clust.profile[[ref.col]])
    }
    clust.nf <- unlist(clust.nf)
    clust.nf[unlist(indices)] <- clust.nf

    # Returning size factors, rather than normalization factors.
    final.sf <- clust.nf * lib.sizes
    if (rescale) { final.sf <- final.sf/mean(final.sf[final.sf>0]) }
    return(final.sf)
}



