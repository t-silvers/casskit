# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import csv
from dataclasses import dataclass, field
import io
import os
from pathlib import Path
import subprocess
from typing import Optional

import pandas as pd

import casskit.io.base as base
import casskit.io.utils as io_utils
import casskit.config as config


@dataclass
class TCGABiolinksSubtype(base.DataURLMixin):

    cache_dir: Optional[Path] = field(init=True, default=None)
    r_script: str = field(init=False, default=Path(__file__).parent / "tcgabiolinks.R")
    
    @io_utils.cache_on_disk
    def fetch(self):
        return self.query_tcgabiolinks()

    def query_tcgabiolinks(self):
        subtypes_l = []
        with subprocess.Popen(["Rscript", self.r_script], stdout=subprocess.PIPE) as p:
            with io.TextIOWrapper(p.stdout, newline=os.linesep) as f:
                reader = csv.reader(f, delimiter=",")
                for r in reader:
                    subtypes_l.append(pd.Series(r))
        
        return pd.concat(subtypes_l[2:], axis=1).set_index(0).T

    def set_cache(self, cache_dir):
        self.path_cache = Path(cache_dir, f"tcga_subtype_tcgabiolinks.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @classmethod
    def get_data(cls, cache_only: bool = False):
        data = cls().fetch()
        if cache_only is False:
            return data
        
    def __post_init__(self):
        if self.cache_dir is None:
            self.cache_dir = Path(config.CACHE_DIR)

        self.set_cache(self.cache_dir)
        self.r_script = Path(self.r_script).as_posix()
        
get_subtypes = TCGABiolinksSubtype.get_data
"""Convenience function for TCGA subtypes from TCGAbiolinks."""
