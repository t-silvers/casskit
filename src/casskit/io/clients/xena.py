# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd
import requests


class XenaClient: ...

class PCAWGXenaLoader(XenaClient):
    def __init__(
        self,
        pcawg_data: PCAWGData,
        cache_dir: Optional[Path] = CACHE_DIR,
    ):
        check_package_version("pyarrow")
        self.omic = pcawg_data.omic
        self.stem = pcawg_data.stem
        self.set_cache(cache_dir)
        self.raw_data = self.fetch()

    @property
    def url(self):
        return f"https://pcawg-hub.s3.us-east-1.amazonaws.com/download/{self.stem}"

    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        return self._fetch(self.url)

    @DataURLMixin.safe_fetch
    def _fetch(self, url) -> pd.DataFrame:
        print("Fetching from, ", url)
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
        cache_dir: Path = CACHE_DIR,
        overwrite: bool = False
    ) -> None:
        for pcawg_data in PCAWG_XENA_DATASETS.values():
            print(f"Building {pcawg_data}")
            cls(pcawg_data, cache_dir).raw_data

    def __repr__(self):
        return "XenaClient"