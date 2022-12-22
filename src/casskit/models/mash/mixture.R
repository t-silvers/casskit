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
    make_option("--random", type = "character",
      metavar = "character"),
    make_option(c("-c", "--covariance"), type = "character",
      help = "Data-driven covariance file", metavar = "character"),
    make_option(c("-o", "--out"), type = "character",
      help = "Output file", metavar = "character")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# Read in data
data.random <- readr::read_rds(opt$random)
U.ed <- readr::read_rds(opt$covariance)

print('Fitting mash model (estimate mixture proportions) ...')
U.c <- mashr::cov_canonical(data.random)
mash.random <- mashr::mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)
mix.proportions <- ashr::get_fitted_g(mash.random)

# Save RDS
readr::write_rds(mix.proportions, opt$out)