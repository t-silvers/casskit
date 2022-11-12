from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

import casskit.utils as utils
import casskit.config as config
import casskit.io.base as base
import casskit.io.utils as io_utils


@dataclass
class BioGRID(base.DataURLMixin):
    """Fetch BioGRID interaction data."""

    cache_dir: Optional[Path] = field(init=True, default=config.get_cache())
    url: str = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.212/BIOGRID-ALL-4.4.212.tab3.zip"
    organism: str = "Homo sapiens"
    
    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:
        return (pd.read_csv(self.url, compression="zip", sep="\t")
                .pipe(utils.column_janitor)
                .query("""
                       organism_name_interactor_a == @self.organism and \
                       organism_name_interactor_b == @self.organism
                       """)
                .filter([
                    "official_symbol_interactor_a", "official_symbol_interactor_b",
                    "synonyms_interactor_a", "synonyms_interactor_b",
                ]))

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"biogrid.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @classmethod
    def get(cls) -> pd.DataFrame:
        return cls().fetch()

    def __post_init__(self):
        self.set_cache(self.cache_dir)

get_biogrid = BioGRID.get
"""Convencience functions for loading BioGRID data."""