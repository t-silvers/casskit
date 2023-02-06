# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from casskit.io.pcawg.config import *
from casskit.io import base
import casskit.io.utils as io_utils
from casskit import descriptors
from casskit import config


class PCAWGXenaLoader(base.DataURLMixin):
    pcawg_data = descriptors.OneOf(*PCAWG_XENA_DATASETS)

    def __init__(
        self,
        pcawg_data: PCAWGData,
        cache_dir: Optional[Path] = None,
    ):
        io_utils.check_package_version("pyarrow")
        
        self.omic = pcawg_data.omic
        self.stem = pcawg_data.stem

        if cache_dir is None:
            cache_dir = config.get_cache()

        self.set_cache(cache_dir)
        self.raw_data = self.fetch()

    @property
    def url(self):
        return f"https://pcawg-hub.s3.us-east-1.amazonaws.com/download/{self.stem}"

    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:        
        return self._fetch(self.url)

    @base.DataURLMixin.safe_fetch
    def _fetch(self, url) -> pd.DataFrame:
        return pd.read_csv(url, sep="\t")
    
    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"PCAWG.{self.omic}.raw.parquet")
        self.read_cache = lambda cache: pd.read_parquet(cache)
        self.write_cache = lambda data, cache: data.to_parquet(cache, engine="pyarrow")

    @classmethod
    def get(cls, data: str) -> pd.DataFrame:
        return cls(PCAWG_XENA_DATASETS[data]).raw_data

    @classmethod
    def build_cache(
        cls,
        cache_dir: Path = None,
        overwrite: bool = False
    ) -> None:
        for pcawg_data in PCAWG_XENA_DATASETS.values():
            print(f"Building {pcawg_data}")
            if cache_dir is None:
                cache_dir = config.get_cache()

            cls(pcawg_data, cache_dir).raw_data

get_pcawg = PCAWGXenaLoader.get
"""Shortcut for TCGAXenaLoader.build_cache"""

build_pcawg = PCAWGXenaLoader.build_cache
"""Shortcut for TCGAXenaLoader.build_cache"""

class PCAWGRawData:
    def __init__(self):
        for data_set in [
            "copynumber",
            "rnaseq",
            "phenotype"
        ]:
            setattr(self, data_set, get_pcawg(data_set))
    
    def __repr__(self):
        return "PCAWGRawData"

class PCAWGDataSet:
    def __init__(self, ret_union: bool = False):
        self.raw_data = PCAWGRawData()
        
        # Prepare data
        copynumber = self.configure_copynumber(self.raw_data.copynumber)
        expression = self.configure_expression(self.raw_data.rnaseq)
        phenotype = self.configure_phenotype(self.raw_data.phenotype)        
        
        if ret_union is True:
            copynumber_sample_ids = copynumber["sample"].unique()
            expression_sample_ids = expression.columns
            shared_samples = np.intersect1d(copynumber_sample_ids,
                                            expression_sample_ids,
                                            assume_unique=True)

            copynumber = copynumber.query("sample in @shared_samples")
            expression = expression.loc[:, shared_samples]
            phenotype = phenotype.query("sample in @shared_samples")
            
        self.copynumber = copynumber
        self.expression = expression
        self.phenotype = phenotype

    @staticmethod
    def configure_copynumber(data):
        return data.rename(columns={"sampleID": "sample"})

    @staticmethod
    def configure_expression(data):
        return (data
                # Parse gene IDs to remove periods
                .assign(gene_id=lambda x: \
                    x.feature.str.split(".", expand=True).get(0))
                .set_index("gene_id")
                .drop("feature", axis=1))

    @staticmethod
    def configure_phenotype(data):
        return (data
                .assign(cancer=lambda df: \
                    io_utils.translate_pcawg_cancer(df.dcc_project_code))
                .explode("cancer")
                .dropna()
                .rename(columns={"icgc_specimen_id": "sample"}))

    def __repr__(self):
        return "PCAWGDataSet"
