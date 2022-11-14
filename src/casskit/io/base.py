from __future__ import annotations

from abc import ABC, abstractmethod
from contextlib import suppress
from functools import wraps
import logging
from pathlib import Path
from typing import Callable, Optional
import urllib
from urllib.request import Request, urlopen

import pandas as pd

import casskit.io.utils as io_utils


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

    def safe_fetch(f: Callable) -> Callable:
        @wraps(f)
        def wrapper(self, *args, **kwargs):
            with suppress(urllib.error.HTTPError):
                return f(self, *args, **kwargs)
        return wrapper

class ElsevierLink(DataURLMixin):
    """Data parser for Springer, Elsevier."""
    def __init__(
        self,
        url,
        skiprows: Optional[int] = 1,
        cache_name: Optional[str] = None,
        cache_dir: Optional[Path] = Path("~/.cache").expanduser(),
    ):
        io_utils.check_package_version("openpyxl")
        self.url = url
        self.skiprows = skiprows        
        self.set_cache(cache_dir, cache_name)

    @property
    def raw_data(self):
        return self.fetch()

    @io_utils.cache_on_disk
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
