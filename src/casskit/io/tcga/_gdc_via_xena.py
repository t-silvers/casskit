from __future__ import annotations

from collections import namedtuple
import logging
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from casskit.io import DataURLMixin
from casskit.io._utils import cache_on_disk
from casskit.descriptors import OneOf
from ...config import CACHE_DIR


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

TCGA_XENA_DATASETS = {
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


class TCGAXenaLoader(DataURLMixin):
    
    cancer = OneOf(*TCGA_CANCERS)
    xena_data = OneOf(*TCGA_XENA_DATASETS)

    def __init__(
        self,
        cancer: str,
        xena_data: XenaData,
        cache_dir: Optional[Path] = CACHE_DIR,
    ):
        self.cancer = cancer
        self.omic = xena_data.omic
        self.sep = xena_data.sep
        self.compression = xena_data.compression

        self.stem = f"{self.cancer}.{self.omic}"
        self.set_cache(cache_dir)
        self.raw_data = self.fetch()
        self.units = self.metadata.loc["unit"][0] if xena_data.units == "from_metadata" else None

    @property
    def basename(self):
        return f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/{self.stem}"

    @property
    def metadata(self) -> Dict[str, object]:
        return pd.read_json(f"{self.basename}.tsv.json", orient="index")

    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        url = self.basename + ".tsv"
        if self.compression == "gzip":
            url += ".gz"
        try:
            return pd.read_csv(url, sep=self.sep, compression=self.compression)
        except Exception as e:
            logging.error(f"Error fetching {url}: {e}")
            raise
    
    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"{self.stem}.raw.parquet")
        self.read_cache = lambda cache: pd.read_parquet(cache)
        self.write_cache = lambda data, cache: data.to_parquet(cache, engine="pyarrow")

    @classmethod
    def get(cls, cancer: str, data: str, cache_dir: Path) -> pd.DataFrame:
        return cls(cancer, TCGA_XENA_DATASETS[data], cache_dir).raw_data

    @classmethod
    def build_cache(cls, cancer: str, cache_dir: Path, overwrite: bool = False) -> None:
        for xena_data in TCGA_XENA_DATASETS.values():
            print(f"Building {xena_data}")

            # exceptions to patterns
            gzip_exceptions = (xena_data.omic == "GDC_phenotype" or cancer == "GDC-PANCAN")
            compression = xena_data.compression if not gzip_exceptions else "gzip"

            # Make new tuple with compression
            xd_ = XenaData(xena_data.omic, xena_data.sep, compression, xena_data.units)
            cls(cancer, xd_, cache_dir).raw_data

get_tcga = TCGAXenaLoader.get
"""Shortcut for TCGAXenaLoader.build_cache"""

build_tcga = TCGAXenaLoader.build_cache
"""Shortcut for TCGAXenaLoader.build_cache"""