from pathlib import Path
from typing import Optional

import pandas as pd

import casskit.io.base as base
import casskit.config as config


class TCGAdataCDRSurvivalLiu2018(base.ElsevierLink):
    """TCGA clinical data resource.
    
    An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Liu et al. 2018.
    DOI: https://doi.org/10.1016/j.cell.2018.02.052
    """
    
    URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx"
    
    def __init__(self, cache_dir: Optional[Path] = None):
        if cache_dir is None:
            cache_dir = Path(config.CACHE_DIR)

        super().__init__(
            url=self.URL,
            skiprows=None,
            cache_name="tcga_cdr_survival",
            cache_dir=cache_dir,
        )
        
    @classmethod
    def get_data(cls, **kwargs) -> pd.DataFrame:
        return cls(**kwargs).raw_data

get_tcga_cdr_survival = TCGAdataCDRSurvivalLiu2018.get_data
"""Convenience function for survival data from Liu et al. 2018."""