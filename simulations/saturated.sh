cat brittlesim.R | sed "s/true.means <- 2^runif(ngenes, 3, 6)/true.means <- 2^runif(ngenes, 9, 12)/" |
                sed 's/output.dir <- "results"/output.dir <- "results_high"/' > highcounts.R
Rscript highcounts.R
