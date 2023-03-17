# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..utils import cache_on_disk


# TODO: Replace with pyrty function
SUBTYPES_ASSET = pd.read_csv(Path(__file__).parent / "assets/subtypes.tsv", sep="\t")

@dataclass
class TCGABiolinksSubtype(DataURLMixin):

    cache_dir: Optional[Path] = field(init=True, default=None)
    
    @cache_on_disk
    def fetch(self):
        return SUBTYPES_ASSET

    def set_cache(self, cache_dir):
        self.path_cache = Path(cache_dir, f"tcga_subtype_tcgabiolinks.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @classmethod
    def get_data(cls, cache_only: bool = False, data: str = None):
        data = cls().fetch()
        if cache_only is False:
            return data
        
    def __post_init__(self):
        if self.cache_dir is None:
            self.cache_dir = CACHE_DIR

        self.set_cache(self.cache_dir)

get_subtypes = TCGABiolinksSubtype.get_data
"""Convenience function for TCGA subtypes from TCGAbiolinks."""
