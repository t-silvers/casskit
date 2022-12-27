#!/usr/bin/env Rscript
# Path: cneqtl/src/cneqtl/mash/summary.R
# Author: Thomas R Silvers <thomas.silvers.1@gmail.com>
# Date: 2022-12-15
# Description: Get posterior estimates and BFs from mash model
# Usage: 

library(dplyr)
library(magrittr)
library(mashr)
library(optparse)
library(readr)
library(tibble)

option_list <- list( 
    make_option("--data", type = "character",
      metavar = "character"),
    make_option("--mix", type = "character",
      metavar = "character")
    )

opt <- parse_args(OptionParser(option_list=option_list))


prepare.mash.output <- function(mash.model){

  pm.res <-
    tibble::as_tibble(mash.model$result$PosteriorMean) %>%
    dplyr::rename_with(~paste0(., "_pm"))

  psd.res <-
    tibble::as_tibble(mash.model$result$PosteriorSD) %>%
    dplyr::rename_with(~paste0(., "_psd"))

  lfsr.res <-
    tibble::as_tibble(mash.model$result$lfsr) %>%
      dplyr::rename_with(~paste0(., "_lfsr"))

  bf.res <-
    tibble::as_tibble(mashr::get_log10bf(mash.model)[,1]) %>%
    dplyr::rename("log10bf"=value)

  # Combine results and return
  dplyr::bind_cols(list(pm.res, psd.res, lfsr.res)) %>%
  dplyr::mutate(model_fit_loglik = mash.model$loglik)
}

# Read in data
mash.data <- readr::read_rds(opt$data)
mix.proportions <- readr::read_rds(opt$mix)

print('Computing posterior summary ...')
m.c <- mashr::mash(mash.data, g = mix.proportions, fixg = TRUE)
mash.res <- prepare.mash.output(m.c)

# Write output
try(writeLines(readr::format_csv(mash.res), stdout()), silent=TRUE)
