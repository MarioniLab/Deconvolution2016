cat brittlesim.R | sed "s/true.means <- 2^runif(ngenes, 3, 6)/true.means <- 2^runif(ngenes, 9, 12)/" |
                sed 's/output.dir <- "results_mid"/output.dir <- "results_high"/' > highcounts.R
Rscript highcounts.R

cat brittlesim.R | sed "s/true.means <- 2^runif(ngenes, 3, 6)/true.means <- rgamma(ngenes, 2, 2)/" |
                sed 's/output.dir <- "results_mid"/output.dir <- "results_low"/' | 
                sed 's/dispersions <- .*/dispersions <- 0.1/' > lowcounts.R
Rscript lowcounts.R
                
