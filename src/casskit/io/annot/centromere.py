from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from casskit.io import base
from casskit import config
import casskit.io.utils as io_utils


@dataclass
class Centromere(base.DataURLMixin):
    UCSC_URLS = {
        "hg19": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz",
        "hg37": "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz",
        "hg38": "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz",
    }
    
    assembly: str
    cache_dir: Optional[Path] = None
    data: pd.DataFrame = field(init=False)
    
    def set_cache(self, cache_dir):        
        self.path_cache = Path(cache_dir, f"centromeres_{self.assembly}.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @io_utils.cache_on_disk
    def fetch(self):

        url = self.UCSC_URLS[self.assembly]
        usecols = [1, 2, 3]
        names = ["chrom", "chromStart", "chromEnd"]
        
        if self.assembly in ["hg19", "hg37"]:
            usecols += [7]
            names += ["gap_type"]        
        
        cen_data = self._fetch(url, usecols=usecols, names=names)
        
        if self.assembly == "hg38":
            return (cen_data
                    .groupby('chrom', as_index=False)
                    .agg({'chromStart': 'min', 'chromEnd': 'max'}))

        elif self.assembly in ["hg19", "hg37"]:
            return (cen_data
                    .query("gap_type == 'centromere'")
                    .drop("gap_type", axis=1)
                    .reset_index(drop=True))
    
    @staticmethod
    def _fetch(url, usecols, names):
        return pd.read_csv(url, sep='\t', usecols=usecols,
                           names=names, compression="gzip")

    def __post_init__(self):
        if self.cache_dir is None:
            self.cache_dir = Path(config.CACHE_DIR)

        self.set_cache(self.cache_dir)
        
        # Parse data
        self.data = self.fetch()