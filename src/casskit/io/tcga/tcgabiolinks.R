pkgs <- c("readr", "TCGAbiolinks")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)

subtypes <- PanCancerAtlas_subtypes()
try(writeLines(readr::format_csv(subtypes), stdout()), silent=TRUE)
