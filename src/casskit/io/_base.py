from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional
from urllib.request import Request, urlopen

import pandas as pd

from _utils import cache_on_disk, check_package_version


class DataURLMixin(ABC):
    """Mixin for data loaders that fetch data from a URL."""

    @property
    def header(self):
        return {
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_14_6) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/14.1.2 Safari/605.1.15",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8"
        }
    
    @abstractmethod
    def fetch(self) -> pd.DataFrame:
        pass

    @abstractmethod
    def set_cache(self):
        pass

class ElsevierLink(DataURLMixin):
    """Data parser for Springer, Elsevier."""
    def __init__(
        self,
        url,
        skiprows: Optional[int] = 1,
        cache_name: Optional[str] = None,
        cache_dir: Optional[Path] = Path("~/.cache").expanduser(),
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
