# Using deconvolution to normalize scRNA-seq data with many zeroes

## Overview

This repository contains the code and manuscript files for the paper **Pooling across cells to normalize single-cell RNA sequencing data with many zero counts**,
by [Lun et al. (2016)](https://doi.org/10.1186/s13059-016-0947-7).

**Note:** Further updates and development of the analysis and simulation code will take place at https://github.com/MarioniLab/FurtherNorm2018. 
If you have general questions regarding the code (i.e., not specifically involving the manuscript), please post your issues at the above repository instead.

## Simulations

To run the simulation code, enter `simulations` and then:

1. Run `lowcounts.R` to perform the low-count simulations, or `brittlesim.R` to perform the high-count simulations.
2. Run `standerr.R` to estimate the variance of the size factor estimates across methods.
3. Run `poolsim.R` to compare the variability of the estimates with and without the ring arrangement.
4. Run `complexity.R` to determine the time-complexity of the deconvolution method.

You can also run `fewcounts.R` to see behaviour with few cells, or `highcounts.R` to see behaviour at very high counts.

## Real analyses

To run the real data analysis code:

1. Make a `data` subdirectory and download the Zeisel et al. tables (http://linnarssonlab.org/cortex) and the Klein data (supplementary tables in GSM1599494, GSM1599499).
2. Enter the `realdata` directory and run `Zeisel.R` and `Klein.R` to pre-process the data and estimate size factors for all cells in each of those two data sets.
3. Run `edgeR.R` to identify DE genes in each data set, and `GOAnalysis.R` to perform a GO analysis on the DE genes.
4. Run `HVGAnalysis.R` to identify highly variable genes in each data set.
5. Run `switchTestedgeR.R` to perform the offset/covariate switching analysis.

Also, run `plotKleinParam.R` to generate plots that justify parameter settings in the simulations.

## Manuscript

The `manuscript` directory contains all LaTeX code used to generate the manuscript.
This can be compiled with `make`.
It assumes that all of the simulations and real data analyses have already been performed.
