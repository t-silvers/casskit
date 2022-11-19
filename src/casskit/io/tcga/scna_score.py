from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


@dataclass
class SCNAScoreXian2021:
    """Single-sample SCNA scores.
    
    "A single somatic copy number alteration (SCNA) score inclusive of
    whole-chromosome, chromosome arm, and focal alterations in a pan-cancer
    analysis of 9,375 samples in The Cancer Genome Atlas (TCGA) database"
    https://www.embopress.org/doi/full/10.15252/embr.202152509
    """
    
    url = "https://github.com/cartercompbio/SCNA_score_analysis/raw/master/data/SCNA/SCNADf_v3.tsv"
    data: pd.DataFrame = field(init=False, default=None)
    
    @classmethod
    def load(cls):
        return cls().data
    
    def __post_init__(self):
        self.data = pd.read_csv(self.url, sep="\t", index_col=0)
        
get_scna_scores = SCNAScoreXian2021.load()