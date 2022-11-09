pkgs <- c("optparse", "readr", "TCGAbiolinks")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE, quietly = TRUE)

# opt_list = list(
#   make_option(c("-d", "--cachedir"), default = getwd(), type = "character",
#               help = "Directory to store downloaded files")
# )

# opt = optparse::parse_args(optparse::OptionParser(option_list=opt_list))

subtypes <- PanCancerAtlas_subtypes()
# readr::write_tsv(subtypes, file.path(opt$cachedir, "subtypes.tsv"))
subtypes