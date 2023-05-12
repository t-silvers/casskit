from __future__ import annotations

from pathlib import Path
from typing import Optional
from urllib.request import Request, urlopen

import pandas as pd

from .base import DataURLMixin
from .config import CACHE_DIR
from .utils import cache_on_disk, check_package_version


class ElsevierLink(DataURLMixin):
    """Data parser for Springer, Elsevier."""
    def __init__(
        self,
        url,
        skiprows: Optional[int] = 1,
        cache_name: Optional[str] = None,
        cache_dir: Optional[Path] = CACHE_DIR,
    ):
        check_package_version("openpyxl")
        self.url = url
        self.skiprows = skiprows        
        self.set_cache(cache_dir, cache_name)

    @property
    def raw_data(self):
        return self.fetch()

    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        request = Request(self.url, None, headers=self.header)
        with urlopen(request) as response:
            return pd.read_excel(response.read(),
                                 engine="openpyxl",
                                 skiprows=self.skiprows)

    def set_cache(self, cache_dir: Path, cache_name: str = None) -> Path:
        if cache_name is None:
            print(self)

        self.path_cache = Path(cache_dir, f"{cache_name}.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)
