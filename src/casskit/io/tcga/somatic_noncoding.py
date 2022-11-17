# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import certifi
from pathlib import Path
import ssl
from typing import Optional
from urllib.request import Request, urlopen

import pandas as pd

import casskit.io.base as base
import casskit.io.utils as io_utils
import casskit.config as config


class SomaticNoncodingZhang2018(base.DataURLMixin):
    """Somatic noncoding mutations from WGS
    
    "Somatic noncoding mutations from 930 tumors were called as described in the main text."
    
    """
    NCMUTS_URL = "https://idekerlab.ucsd.edu/papers/wzhang2017/Somatic_mutations_TCGA_930.txt"
    COVARS_URL = "https://idekerlab.ucsd.edu/papers/wzhang2017/covariates.txt"
    
    def __init__(
        self,
        cache_dir: Optional[Path] = None,
    ):
        if cache_dir is None:
            cache_dir = Path(config.CACHE_DIR)

    @property
    def prepared_data(self) -> pd.DataFrame:
        return self.prepare()
    
    @io_utils.cache_on_disk
    def fetch(self) -> pd.DataFrame:
        
        # This hangs
        context = ssl.create_default_context(cafile=certifi.where())
        context.check_hostname = False
        context.verify_mode = ssl.CERT_NONE

        request = Request(self.NCMUTS_URL, None)
        with urlopen(request, context=context) as response:
            return pd.read_table(response.read())

    @classmethod
    def get_data(cls, **kwargs) -> pd.DataFrame:
        return cls(**kwargs).prepared_data

get_zhang_noncoding = SomaticNoncodingZhang2018.get_data
"""Convenience function for somatic noncoding mutations from Zhang et al. 2018."""
