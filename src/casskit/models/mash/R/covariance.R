#!/usr/bin/env Rscript
# Path: cneqtl/src/cneqtl/mash/covariance.R
# Author: Thomas R Silvers <thomas.silvers.1@gmail.com>
# Date: 2022-12-15
# Description: Calculate mash data-driven covariances
# Usage: 
# Notes: https://stephenslab.github.io/mashr/articles/eQTL_outline.html

library(dplyr)
library(magrittr)
library(mashr)
library(optparse)
library(readr)
library(tibble)

option_list <- list( 
    make_option("--strong", type = "character",
      metavar = "character"),
    make_option(c("-o", "--out"), type = "character",
      help = "Output file", metavar = "character")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# Read in data
data.strong <- readr::read_rds(opt$strong)

print('Using strong tests to set up data-driven covariances ...')
npcs.to.use <- 5 # number of PCs used in https://stephenslab.github.io/mashr/articles/eQTL_outline.html
if (dim(data.strong$Bhat)[2] <= npcs.to.use) {
  # Error: npc not less than or equal to n_conditions(data)
  npcs.to.use <- dim(data.strong$Bhat)[2] - 1
}
U.pca <- mashr::cov_pca(data.strong, npcs.to.use)
U.ed <- mashr::cov_ed(data.strong, U.pca)

# Save RDS file
readr::write_rds(U.ed, opt$out)