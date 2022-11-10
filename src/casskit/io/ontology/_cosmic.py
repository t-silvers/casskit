from dataclasses import dataclass
import logging
from pathlib import Path
import subprocess
from this import d

import pandas as pd

from .._utils import cache_on_disk


@dataclass
class CancerGeneCensusCOSMIC:
    """COSMIC Cancer Gene Census.
    
    https://cancer.sanger.ac.uk/cosmic/download
    
    Notes
    -----
    Users must have a registered account with COSMIC to download. Please use
    your own account credentials to access.
    """
    
    cache_dir: Path = Path("~/.cache").expanduser()
    acct_email: str = "tsilvers@stanford.edu"
    acct_pwd: str = "Public212Password&"
    
    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        cmd = (
            '''
            AUTHENTICATION=$(echo "email:pwd" | base64)
            LINK=$(curl -H "Authorization: Basic ${AUTHENTICATION}" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cosmic/v96/cancer_gene_census.csv | jq -r '.url')
            curl $LINK >> cache
            '''
            .replace("email", self.acct_email)
            .replace("pwd", self.acct_pwd)
            .replace("cache", self.path_cache.as_posix())
        )

        try:
            subprocess.check_output(cmd, shell=True)
            return pd.read_csv(self.path_cache)

        except Exception as e:
            logging.error(f"Error fetching data using {cmd}: {e}")
            raise

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, "cosmic.csv")
        self.read_cache = lambda cache: pd.read_csv(cache)
        self.write_cache = lambda data, cache: data.to_csv(cache, index=False)

    @classmethod
    def get(cls, cache_dir: Path = Path("~/.cache").expanduser()) -> pd.DataFrame:
        return cls(cache_dir).fetch()

    def __post_init__(self):
        self.set_cache(self.cache_dir)

get_cosmic = CancerGeneCensusCOSMIC.get
"""Convencience functions for loading COSMIC data."""