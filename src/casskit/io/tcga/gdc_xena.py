# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from .config import *
import casskit.config as config
from casskit.descriptors import OneOf
from casskit.io.base import DataURLMixin
import casskit.io.utils as io_utils


class TCGAXenaLoader(DataURLMixin):
    
    cancer = OneOf(*TCGA_CANCERS)
    xena_data = OneOf(*tcga_xena_datasets)

    def __init__(
        self,
        cancer: str,
        xena_data: XenaData,
        cache_dir: Optional[Path] = None,
        minimal: bool = False,
        num_samples: int = 50,
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
        
        # For testing and memory-limited settings
        self.minimal = minimal
        self.num_samples = num_samples
        self.raw_data = self.fetch()

    @property
    def basename(self):
        return f"https://gdc-hub.s3.us-east-1.amazonaws.com/download/{self.stem}"

    @property
    @DataURLMixin.safe_fetch
    def metadata(self) -> TCGAXenaMetadata:
        metadata_df = pd.read_json(f"{self.basename}.tsv.json", orient="index")
        return TCGAXenaMetadata(metadata_df)

    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:
        url = self.basename + ".tsv"
        if self.compression == "gzip":
            url += ".gz"
        
        return self._fetch(url)

    @DataURLMixin.safe_fetch
    def _fetch(self, url) -> pd.DataFrame:
        data = pd.read_csv(url, sep=self.sep, compression=self.compression)
        if self.minimal is True:
            
            # Get minimal samples
            if self.metadata.type == "genomicMatrix":
                data = (data
                        .set_index(data.columns[0])
                        .sort_index(axis=1)
                        .iloc[:, :self.num_samples]
                        .reset_index())

            elif self.metadata.type == "genomicSegment":
                samples = data["sample"].unique().tolist()
                samples = sorted(samples)[:self.num_samples]
                data = data.query("sample in @samples")

            elif ((self.metadata.type == "clinicalMatrix") or
                  (self.metadata.type == "mutationVector")):
                samples = data.iloc[:, 0].unique().tolist()
                samples = sorted(samples)[:self.num_samples]
                data = data[data[data.columns[0]].isin(samples)]

        return data
    
    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"GDC.{self.cancer}/data/{self.omic}.raw.parquet")
        self.read_cache = lambda cache: pd.read_parquet(cache)
        self.write_cache = lambda data, cache: data.to_parquet(cache, engine="pyarrow")

    @classmethod
    def get(cls, cancer: str, data: str) -> pd.DataFrame:
        return cls(cancer, tcga_xena_datasets[data]).raw_data

    @classmethod
    def build_cache(
        cls,
        cancer: str,
        cache_dir: Path = None,
        overwrite: bool = False,
        minimal: bool = False,
    ) -> None:
        for xena_data in tcga_xena_datasets.values():
            print(f"Building {xena_data}")
            if cache_dir is None:
                cache_dir = config.get_cache()

            # exceptions to patterns
            gzip_exceptions = (xena_data.omic == "GDC_phenotype" or cancer == "GDC-PANCAN")
            compression = xena_data.compression if not gzip_exceptions else "gzip"

            # Make new tuple with compression
            xd_ = XenaData(xena_data.omic, xena_data.sep, compression, xena_data.units)
            cls(cancer, xd_, cache_dir=cache_dir, minimal=minimal).raw_data

class TCGAXenaMetadata:
    def __init__(self, data):
        self.data = data
        if self.data is not None:
            for metadata_attr in self.data.index:
                setattr(self, metadata_attr, self.data.loc[metadata_attr][0])

get_gdc_tcga = TCGAXenaLoader.get
"""Shortcut for TCGAXenaLoader.build_cache"""

build_tcga = TCGAXenaLoader.build_cache
"""Shortcut for TCGAXenaLoader.build_cache"""