#!/usr/bin/env Rscript
# Path: cneqtl/src/cneqtl/mash/nullcorr.R
# Author: Thomas R Silvers <thomas.silvers.1@gmail.com>
# Date: 2022-12-15
# Description: Estimate correlation structure across tissues from null tests
# Usage: 
# Notes: https://stephenslab.github.io/mashr/articles/intro_correlations.html
# https://stephenslab.github.io/mashr/reference/estimate_null_correlation_simple.html

library(dplyr)
library(magrittr)
library(mashr)
library(optparse)
library(readr)
library(tibble)

option_list <- list( 
    make_option("--b_rand", type = "character",
      metavar = "character"),
    make_option("--s_rand", type = "character",
      metavar = "character"),
    make_option("--method", type = "character",
      default = "simple", metavar = "character"),
    make_option("--out", type = "character",
      metavar = "character")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# Read in data
bhat.random <- as.matrix(read_csv(opt$b_rand))
shat.random <- as.matrix(read_csv(opt$s_rand))
data.random <- mashr::mash_set_data(bhat.random, shat.random) # mash data without variance

print('Estimating correlation structure in the null from random data ...')
V.simple <- mashr::estimate_null_correlation_simple(data.random)

# For tissues with few mapped cn-eQTLs, random subsets are too small.
# Rather than adjusting sampling, just replace corresponding values in
# rank matrix with 0
V.simple[is.na(V.simple)] <- 0

# Alternative: Reduce z_thresh 2 --> 1 --> ?

if (opt$method == "simple"){ 
  
  # Save RDS file
  readr::write_rds(V.simple, opt$out)

} else if (opt$method == "em"){
  
  data.Vsimple = mashr::mash_update_data(data.random, V=V.simple)
  U.c = mashr::cov_canonical(data.Vsimple) 
  
  V.em = mashr::mash_estimate_corr_em(data.random, U.c, details = TRUE)
  m.Vem = V.em$mash.model

  # Save RDS file
  readr::write_rds(V.em, opt$out)

}