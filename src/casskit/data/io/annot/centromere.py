from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..utils import cache_on_disk


@dataclass
class Centromere(DataURLMixin):
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

    @cache_on_disk
    def fetch(self):
        url = self.UCSC_URLS[self.assembly]
        usecols = [1, 2, 3]
        names = ["Chromosome", "Start", "End"]
        
        if self.assembly in ["hg19", "hg37"]:
            usecols += [7]
            names += ["gap_type"]        
        
        cen_data = self._fetch(url, usecols=usecols, names=names)
        
        if self.assembly == "hg38":
            return (cen_data
                    .groupby("Chromosome", as_index=False)
                    .agg({"Start": "min", "End": "max"}))

        elif self.assembly in ["hg19", "hg37"]:
            return (cen_data
                    .query("gap_type == 'centromere'")
                    .drop("gap_type", axis=1)
                    .reset_index(drop=True))
    
    @staticmethod
    def _fetch(url, usecols, names):
        print(url)
        return pd.read_csv(url, sep='\t', usecols=usecols,
                           names=names, compression="gzip")

    def __post_init__(self):
        if self.cache_dir is None:
            self.cache_dir = CACHE_DIR

        self.set_cache(self.cache_dir)
        
        # Parse data
        self.data = self.fetch()

def annotate_chrom_arm(data, assembly):
    """Annotate chromosome arm."""
    # Input checks
    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame")
    
    if not "Chromosome" in data.columns:
        raise ValueError("data must have a 'Chromosome' column")

    if not "Start" in data.columns:
        raise ValueError("data must have a 'Start' column")

    if not "End" in data.columns:
        raise ValueError("data must have a 'End' column")
    
    # Fetch centromere data
    cen_data = (Centromere(assembly).data
                # Downcast to save memory
                .assign(Chromosome=lambda x: x.Chromosome.astype("category"),
                        cen_start=lambda x: x.pop("Start").astype("float"),
                        cen_end=lambda x: x.pop("End").astype("float")))
    
    # Add arm annotation
    return (data
            .merge(cen_data)
            .assign(arm=lambda df:
                df.apply(lambda x: \
                    ["p"] if (x.Start <= x.cen_start) & (x.End < x.cen_end) else 
                    ["q"] if (x.End >= x.cen_end) & (x.Start > x.cen_start) else
                    ["cen"] if (x.Start > x.cen_start) & (x.End < x.cen_end) else
                    # Spans both arms
                    ["p", "q"], axis=1)
            )
            .explode("arm")
            .assign(arm=lambda x: x.arm.astype("category"))
            .drop(["cen_start", "cen_end"], axis=1))