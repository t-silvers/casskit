# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from .tcga_config import TCGA_CANCERS, TCGA_XENA_DATASETS, XenaData
from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..descriptors import OneOf
from ..utils import cache_on_disk, check_package_version


class TCGAXenaLoader(DataURLMixin):
    
    cancer = OneOf(*TCGA_CANCERS)
    xena_data = OneOf(*TCGA_XENA_DATASETS)

    def __init__(
        self,
        cancer: str,
        xena_data: XenaData,
        cache_dir: Optional[Path] = CACHE_DIR,
        minimal: bool = False,
        n_samples: int = 50,
    ):
        check_package_version("pyarrow")
        
        self.cancer = cancer
        self._units = xena_data.units
        self.omic = xena_data.omic
        self.sep = xena_data.sep
        self.compression = xena_data.compression
        self.stem = f"{self.cancer}.{self.omic}"
        
        # Configure caching
        self.set_cache(cache_dir)
        
        # For testing and memory-limited settings
        self.minimal = minimal
        self.n_samples = n_samples
        self.raw_data = self.fetch()

    @property
    def basename(self):
        return (f"https://gdc-hub.s3.us-east-1.amazonaws.com"
                f"/download/{self.stem}")

    @property
    @DataURLMixin.safe_fetch
    def metadata(self) -> TCGAXenaMetadata:
        metadata_df = pd.read_json(f"{self.basename}.tsv.json", orient="index")
        return TCGAXenaMetadata(metadata_df)

    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        url = self.basename + ".tsv"
        if self.compression == "gzip":
            url += ".gz"
        
        return self._fetch(url)

    @DataURLMixin.safe_fetch
    def _fetch(self, url) -> pd.DataFrame:
        print("Fetching data from ", url)
        data = pd.read_csv(url, sep=self.sep, compression=self.compression)
        if self.minimal is True:
            data = subset_tcga_for_testing(data, self.metadata.type,
                                           self.n_samples)

        return data
    
    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"GDC.{self.cancer}/data/{self.omic}.raw.parquet")
        self.read_cache = lambda cache: pd.read_parquet(cache)
        self.write_cache = lambda data, cache: data.to_parquet(cache, engine="pyarrow")

    @classmethod
    def get(cls, cancer: str, data: str) -> pd.DataFrame:
        return cls(cancer, TCGA_XENA_DATASETS[data]).raw_data

    @classmethod
    def build_cache(
        cls,
        cancer: str,
        cache_dir: Path = CACHE_DIR,
        minimal: bool = False,
        n_samples: int = 20,
    ) -> None:
        for xena_data in TCGA_XENA_DATASETS.values():
            print(f"Building {xena_data}")
            # Exceptions to patterns
            gzip_exceptions = (xena_data.omic == "GDC_phenotype" or cancer == "GDC-PANCAN")
            compression = xena_data.compression if not gzip_exceptions else "gzip"

            # Make new tuple with correct compression
            xd_ = XenaData(xena_data.omic, xena_data.sep,
                           compression, xena_data.units)
            
            cls(cancer, xd_, cache_dir=cache_dir,
                minimal=minimal, n_samples=n_samples
                ).raw_data

class TCGAXenaMetadata:
    def __init__(self, data):
        self.data = data
        if self.data is not None:
            for metadata_attr in self.data.index:
                setattr(self, metadata_attr, self.data.loc[metadata_attr][0])

def subset_tcga_for_testing(data, data_type, n_samples):
    # TODO: Check n_samples is possible
    if data_type == "genomicMatrix":
        return (data
                .set_index(data.columns[0])
                .sort_index(axis=1)
                .iloc[:, :n_samples]
                .reset_index())

    elif data_type == "genomicSegment":
        samples = data["sample"].unique().tolist()
        samples = sorted(samples)[:n_samples]
        return data.query("sample in @samples")

    elif ((data_type == "clinicalMatrix") or
            (data_type == "mutationVector")):
        samples = data.iloc[:, 0].unique().tolist()
        samples = sorted(samples)[:n_samples]
        return data[data[data.columns[0]].isin(samples)]

get_gdc_tcga = TCGAXenaLoader.get
"""Shortcut for TCGAXenaLoader.build_cache"""

build_tcga = TCGAXenaLoader.build_cache
"""Shortcut for TCGAXenaLoader.build_cache"""
