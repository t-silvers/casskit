from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

import casskit.utils as utils
import casskit.config as config
import casskit.io.base as base
import casskit.io.utils as io_utils


@dataclass
class CORUM(base.DataURLMixin):
    """Fetch CORUM interaction data.
    
    Note
    ----
    - complex names are not unique, but complex IDs are unique
    """

    cache_dir: Optional[Path] = field(init=True, default=None)
    url: str = "https://mips.helmholtz-muenchen.de/corum/download/releases/current/allComplexes.txt.zip"
    organism: str = "Human"
    
    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:
        return (pd.read_csv(self.url, compression="zip", sep="\t")
                .query("Organism == @self.organism")
                .pipe(utils.column_janitor)
                .assign(
                    subunits_gene_name = lambda x: x['subunits_gene_name'].str.split(';'),
                    complexname = lambda x: x['complexname'].apply(utils.janitor),
                )
                .explode(column='subunits_gene_name'))

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"corum.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @classmethod
    def get(cls, cache_dir: Path = config.CACHE_DIR) -> pd.DataFrame:
        return cls(cache_dir).fetch()

    def __post_init__(self):
        if self.cache_dir is None:
            self.cache_dir = Path(config.CACHE_DIR)
        self.set_cache(self.cache_dir)

get_corum = CORUM.get
"""Convencience functions for loading CORUM data."""