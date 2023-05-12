# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..utils import cache_on_disk, column_janitor


@dataclass
class TRRUST(DataURLMixin):
    """Fetch TRRUST transcription factor data."""
    cache_dir: Optional[Path] = field(init=True, default=CACHE_DIR)
    url: str = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"

    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        COLNAMES_FROM_WEB = ['TranscriptionFactor', 'TargetGene',
                             'ModeOfRegulation', 'article-id_pmid']
        return (pd.read_csv(self.url, sep='\t', names=COLNAMES_FROM_WEB)
                .pipe(column_janitor))

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, "trrust.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @staticmethod
    def parse_reg_direction(s: pd.Series) -> pd.Series:
        return s.str.lower().replace({'repression': -1, 'unknown': 0, 'activation': 1})

    @classmethod
    def get(cls, cache_only: bool = False) -> pd.DataFrame:
        data = cls().fetch()
        if cache_only is False:
            return data

    def __post_init__(self):
        self.set_cache(self.cache_dir)

get_trrust = TRRUST.get
"""Convencience functions for loading TRRUST data."""