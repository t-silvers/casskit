# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from pathlib import Path
from typing import Optional

import pandas as pd

import casskit.io.base as base
import casskit.config as config


class PloidyScoresLOHCiani2022(base.ElsevierLink):
    """Aneuploidy scores from Ciani et al. 2022.
    
    TCGA sample information, including LOH, purity, 6 ploidy measures,
    and clinical info.

    Allele-specific genomic data elucidate the role of somatic gain and
    copy-number neutral loss of heterozygosity in cancer.
    Ciani, Fedrizzi, Prandi et al. 2022.
    DOI: https://doi.org/10.1016/j.cels.2021.10.001
    
    Table S1. Study sample information file, related to Figure 1.
    """

    URL = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471221003835-mmc2.xlsx"
    
    def __init__(self, cache_dir: Optional[Path] = None):
        if cache_dir is None:
            cache_dir = Path(config.CACHE_DIR)

        super().__init__(
            url=self.URL,
            skiprows=1,
            cache_name="aneuploidy_scores_loh",
            cache_dir=cache_dir,
        )
        
    @classmethod
    def get_data(cls, **kwargs) -> pd.DataFrame:
        return cls(**kwargs).raw_data

get_tcga_aneuploidy_scores = PloidyScoresLOHCiani2022.get_data
"""Convenience function for aneuploidy scores from Ciani et al. 2022."""