# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from collections import namedtuple
from pathlib import Path
from typing import Dict, Optional

import pandas as pd

import casskit.io.base as base
import casskit.io.utils as io_utils
import casskit.descriptors as descriptors
import casskit.config as config



PCAWGData = namedtuple("pcawg", ["stem", "omic"])

PCAWG_XENA_DATASETS = {
    "rnaseq": PCAWGData("tophat_star_fpkm_uq.v2_aliquot_gl.sp.log", "rnaseq"),
    "copynumber": PCAWGData("20170119_final_consensus_copynumber_sp", "copynumber"),
    "phenotype": PCAWGData("project_code_sp", "phenotype"),
}


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