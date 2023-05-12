# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..utils import cache_on_disk, column_janitor, janitor


@dataclass
class CORUM(DataURLMixin):
    """Fetch CORUM interaction data.
    
    Note
    ----
    - complex names are not unique, but complex IDs are unique
    """

    cache_dir: Optional[Path] = field(init=True, default=CACHE_DIR)
    url: str = "https://mips.helmholtz-muenchen.de/corum/download/releases/current/allComplexes.txt.zip"
    organism: str = "Human"
    
    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        return (pd.read_csv(self.url, compression="zip", sep="\t")
                .query("Organism == @self.organism")
                .pipe(column_janitor)
                .assign(
                    subunits_gene_name = lambda x: x['subunits_gene_name'].str.split(';'),
                    complexname = lambda x: x['complexname'].apply(janitor),
                )
                .explode(column='subunits_gene_name'))

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, "corum.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @classmethod
    def get(cls, cache_only: bool = False) -> pd.DataFrame:
        data = cls().fetch()
        if cache_only is False:
            return data

    def __post_init__(self):
        self.set_cache(self.cache_dir)

get_corum = CORUM.get
"""Convencience functions for loading CORUM data."""