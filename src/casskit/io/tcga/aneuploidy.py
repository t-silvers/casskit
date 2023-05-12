# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..elsevier import ElsevierLink


CIANI2022_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S2405471221003835-mmc2.xlsx"

@dataclass
class PloidyScoresLOHCiani2022(ElsevierLink):
    """Aneuploidy scores from Ciani et al. 2022.
    
    TCGA sample information, including LOH, purity, 6 ploidy measures,
    and clinical info.

    Allele-specific genomic data elucidate the role of somatic gain and
    copy-number neutral loss of heterozygosity in cancer.
    Ciani, Fedrizzi, Prandi et al. 2022.
    DOI: https://doi.org/10.1016/j.cels.2021.10.001
    
    Table S1. Study sample information file, related to Figure 1.
    
    ascores = ElsevierLink(url="https://ars.els-cdn.com/content/image/1-s2.0-S2405471221003835-mmc2.xlsx", skiprows=1)
    """

    data: pd.DataFrame = field(init=False)
    cache_dir: Optional[Path] = field(init=True, default=None)
    url: str = field(init=False, default=CIANI2022_URL)

    @classmethod
    def get_data(cls, cache_only: bool = False, data: str = None):
        data = cls().data
        if cache_only is False:
            return data

    def __post_init__(self):
        super().__init__(
            self.url,
            skiprows=1,
            cache_name="aneuploidy_scores_loh",
            cache_dir=self.cache_dir
        )
        self.data = self.fetch()

get_tcga_aneuploidy_scores = PloidyScoresLOHCiani2022.get_data
"""Convenience function for aneuploidy scores from Ciani et al. 2022."""