#!/usr/bin/env Rscript
# Path: cneqtl/src/cneqtl/mash/mashdata.R
# Author: Thomas R Silvers <thomas.silvers.1@gmail.com>
# Date: 2022-12-15
# Description: Prepare mash data object
# Usage: 
# Notes: https://stephenslab.github.io/mashr/articles/eQTL_outline.html

library(dplyr)
library(magrittr)
library(mashr)
library(optparse)
library(readr)
library(tibble)

option_list <- list( 
    make_option(c("-b", "--beta"), type = "character",
      metavar = "character"),
    make_option(c("-s", "--sd"), type = "character",
      metavar = "character"),
    make_option("--nullcorr", type = "character",
      help = "Null correlation file", metavar = "character"),
    make_option("--out", type = "character",
      help = "Output strong file", metavar = "character")
    )

opt <- parse_args(OptionParser(option_list=option_list))

# Read in data
bhat <- as.matrix(read_csv(opt$beta))
shat <- as.matrix(read_csv(opt$sd))
Vhat <- readr::read_rds(opt$nullcorr)

print("Setting up main data objects with correlation structure ...")
mash.data <- mashr::mash_set_data(bhat, shat, V=Vhat)
readr::write_rds(mash.data, opt$out)