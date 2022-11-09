pkgs <- c("readr", "TCGAbiolinks")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)

subtypes <- PanCancerAtlas_subtypes()
writeLines(readr::format_csv(subtypes), stdout())