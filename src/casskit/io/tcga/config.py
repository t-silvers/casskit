# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from collections import namedtuple

__all__ = ["TCGA_CANCERS", "XenaData", "tcga_xena_datasets"]


TCGA_CANCERS = [
    "TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-CHOL",
    "TCGA-COAD", "TCGA-DLBC", "TCGA-ESCA", "TCGA-GBM", "TCGA-HNSC",
    "TCGA-KICH", "TCGA-KIRC", "TCGA-KIRP", "TCGA-LAML", "TCGA-LIHC",
    "TCGA-LUAD", "TCGA-LUSC", "TCGA-MESO", "TCGA-OV", "TCGA-PAAD",
    "TCGA-PCPG", "TCGA-PRAD", "TCGA-READ", "TCGA-SARC", "TCGA-SKCM",
    "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-THYM", "TCGA-UCEC",
    "TCGA-UCS", "TCGA-UVM", "GDC-PANCAN"
]

XenaData = namedtuple("xena", ["omic", "sep", "compression", "units"])

tcga_xena_datasets = {
    "cnv": XenaData("cnv", "\t", "gzip", "from_metadata"),
    "GDC_phenotype": XenaData("GDC_phenotype", "\t", "gzip", None),
    "gistic": XenaData("gistic", "\t", "gzip", "from_metadata"),
    "htseq_counts": XenaData("htseq_counts", "\t", "gzip", "from_metadata"),
    "htseq_fpkm": XenaData("htseq_fpkm", "\t", "gzip", "from_metadata"),
    "htseq_fpkm-uq": XenaData("htseq_fpkm-uq", "\t", "gzip", "from_metadata"),
    "masked_cnv": XenaData("masked_cnv", "\t", "gzip", "from_metadata"),
    "methylation27": XenaData("methylation27", "\t", "gzip", "from_metadata"),
    "methylation450": XenaData("methylation450", "\t", "gzip", "from_metadata"),
    "mirna": XenaData("mirna", "\t", "gzip", "from_metadata"),
    "muse_snv": XenaData("muse_snv", "\t", "gzip", None),
    "mutect2_snv": XenaData("mutect2_snv", "\t", "gzip", None),
    "somaticsniper_snv": XenaData("somaticsniper_snv", "\t", "gzip", None),
    "survival": XenaData("survival", "\t", None, None),
    "varscan2_snv": XenaData("varscan2_snv", "\t", "gzip", None),   
}
