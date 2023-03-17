from __future__ import annotations

from abc import ABC, abstractmethod
from contextlib import suppress
from functools import wraps
from typing import Callable
import urllib

import pandas as pd

from .config import HEADER


class DataURLMixin(ABC):
    """Mixin for data loaders that fetch data from a URL."""
    header = HEADER
    
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