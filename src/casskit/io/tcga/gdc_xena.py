# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from collections import namedtuple
import logging
from pathlib import Path
from typing import Dict, Optional
import urllib
import warnings

import pandas as pd

import casskit.io.base as base
import casskit.io.utils as io_utils
import casskit.descriptors as descriptors
import casskit.config as config


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


class TCGAXenaLoader(base.DataURLMixin):
    
    cancer = descriptors.OneOf(*TCGA_CANCERS)
    xena_data = descriptors.OneOf(*TCGA_XENA_DATASETS)

    def __init__(
        self,
        cancer: str,
        xena_data: XenaData,
        cache_dir: Optional[Path] = None,
    ):
        io_utils.check_package_version("pyarrow")
        
        self.cancer = cancer
        self._units = xena_data.units
        self.omic = xena_data.omic
        self.sep = xena_data.sep
        self.compression = xena_data.compression

        self.stem = f"{self.cancer}.{self.omic}"
        if cache_dir is None:
            cache_dir = config.get_cache()

        self.set_cache(cache_dir)
        self.raw_data = self.fetch()

    @property
    def basename(self):
        return f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/{self.stem}"

    @property
    @base.DataURLMixin.safe_fetch
    def metadata(self) -> Dict[str, object]:
        return pd.read_json(f"{self.basename}.tsv.json", orient="index")

    @property
    def units(self):
        if self._units == "from_metadata":
            if (self.metadata is not None & "unit" in self.metadata.index):
                return self.metadata.loc["unit"][0]

        return None

    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:
        url = self.basename + ".tsv"
        if self.compression == "gzip":
            url += ".gz"
        
        return self._fetch(url)

    @base.DataURLMixin.safe_fetch
    def _fetch(self, url) -> pd.DataFrame:
        return pd.read_csv(url, sep=self.sep, compression=self.compression)
    
    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"{self.stem}.raw.parquet")
        self.read_cache = lambda cache: pd.read_parquet(cache)
        self.write_cache = lambda data, cache: data.to_parquet(cache, engine="pyarrow")

    @classmethod
    def get(cls, cancer: str, data: str) -> pd.DataFrame:
        return cls(cancer, TCGA_XENA_DATASETS[data]).raw_data

    @classmethod
    def build_cache(
        cls,
        cancer: str,
        cache_dir: Path = None,
        overwrite: bool = False
    ) -> None:
        for xena_data in TCGA_XENA_DATASETS.values():
            print(f"Building {xena_data}")
            if cache_dir is None:
                cache_dir = config.get_cache()

            # exceptions to patterns
            gzip_exceptions = (xena_data.omic == "GDC_phenotype" or cancer == "GDC-PANCAN")
            compression = xena_data.compression if not gzip_exceptions else "gzip"

            # Make new tuple with compression
            xd_ = XenaData(xena_data.omic, xena_data.sep, compression, xena_data.units)
            cls(cancer, xd_, cache_dir).raw_data

get_gdc_tcga = TCGAXenaLoader.get
"""Shortcut for TCGAXenaLoader.build_cache"""

build_tcga = TCGAXenaLoader.build_cache
"""Shortcut for TCGAXenaLoader.build_cache"""