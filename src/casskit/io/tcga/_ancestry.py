from pathlib import Path
from typing import Optional

import pandas as pd

from casskit.io._base import DataURLMixin
from casskit.io._utils import cache_on_disk

from ...utils import column_janitor
from ...config import CACHE_DIR # TEMP


class TCGAAncestryPCs(DataURLMixin):
    
    SOURCES = {
        "broad": "https://api.gdc.cancer.gov/data/1fee3458-14ee-4b4b-964c-a05164b68066",
        "ucsf": "https://api.gdc.cancer.gov/data/fdfa536a-c3c8-405d-99d9-bc9375b5084c",
        "washu": "https://api.gdc.cancer.gov/data/549b67ce-991d-4356-82fa-d09f9d9a23c8"
    }

    def __init__(
        self,
        cache_dir: Optional[Path] = CACHE_DIR,
        institute: str = "broad",
        sep: str = "\t",
    ):
        self.institute = institute
        self.sep = sep
        self.set_cache(cache_dir)
    
    @cache_on_disk
    def fetch(self):
        return pd.read_csv(self.SOURCES[self.institute], sep=self.sep)

    def set_cache(self, cache_dir):
        self.path_cache = Path(cache_dir, f"tcga_ancestry_{self.institute}.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @staticmethod
    def prepare_ancestry_data(data: pd.DataFrame) -> pd.DataFrame:
        return (data                
                .dropna(how="all", axis=1)
                .set_index("sample")
                .filter(like="pc")
                .reset_index())

    @classmethod
    def get_broad(cls):
        broad_data = cls(institute="broad").fetch()
        broad_data[["pc1", "pc2", "pc3"]] = broad_data.pop("PC1:PC2:PC3").str.split(":", expand=True).astype(float)
        return broad_data.rename(columns={"SampleID": "sample", "Ancestry_assignment": "ancestry"})

    @classmethod
    def get_ucsf(cls):
        return (cls(institute="ucsf", sep=None)
                .fetch()
                .drop(["Unnamed: 0", "ethnicity"], axis=1)
                .rename(columns={
                    "Patient_ID": "sample", "Aliquot_ID": "aliquot_id",
                    "pam.ancestry.cluster": "ancestry", "race": "ethnicity"
                }))

    @classmethod
    def get_washu(cls):
        return (cls(institute="washu")
                .fetch()
                .rename(columns={'Case': 'sample', 'Sample': 'aliquot_id'})
                .pipe(column_janitor))

    @classmethod
    def get_data(cls, cache_only: bool = False):
        broad_data = cls.get_broad()
        ucsf_data = cls.get_ucsf()
        washu_data = cls().get_washu()

        if cache_only is False:
            # PC values change across data set. Use WashU set, as this appears most complete.
            return (cls
                    .prepare_ancestry_data(washu_data)
                    .set_index("sample")
                    .add_prefix("ancestry_"))

get_ancestry_pcs = TCGAAncestryPCs.get_data
"""Convenience function for ancestry PCs from TCGA."""
