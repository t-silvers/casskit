from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

import casskit.utils as utils
import casskit.config as config
import casskit.io.base as base
import casskit.io.utils as io_utils


@dataclass
class TRRUST(base.DataURLMixin):
    """Fetch TRRUST transcription factor data."""

    cache_dir: Optional[Path] = field(init=True, default=None)
    url: str = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"

    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:
        COLNAMES_FROM_WEB = ['TranscriptionFactor', 'TargetGene',
                             'ModeOfRegulation', 'article-id_pmid']
        return (pd.read_csv(self.url, sep='\t', names=COLNAMES_FROM_WEB)
                .pipe(utils.column_janitor))

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, "trrust.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @staticmethod
    def parse_reg_direction(s: pd.Series) -> pd.Series:
        return s.str.lower().replace({'repression': -1, 'unknown': 0, 'activation': 1})

    @classmethod
    def get(cls) -> pd.DataFrame:
        return cls().fetch()

    def __post_init__(self):
        if self.cache_dir is None:
            self.cache_dir = config.CACHE_DIR
        self.set_cache(self.cache_dir)

get_trrust = TRRUST.get
"""Convencience functions for loading TRRUST data."""