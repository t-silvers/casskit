# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import pandas as pd

from .pcawg_config import PCAWG_XENA_DATASETS, PCAWGData
from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..descriptors import OneOf
from ..utils import cache_on_disk, check_package_version


class PCAWGXenaLoader(DataURLMixin):
    pcawg_data = OneOf(*PCAWG_XENA_DATASETS)

    def __init__(
        self,
        pcawg_data: PCAWGData,
        cache_dir: Optional[Path] = CACHE_DIR,
    ):
        check_package_version("pyarrow")
        self.omic = pcawg_data.omic
        self.stem = pcawg_data.stem
        self.set_cache(cache_dir)
        self.raw_data = self.fetch()

    @property
    def url(self):
        return f"https://pcawg-hub.s3.us-east-1.amazonaws.com/download/{self.stem}"

    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        return self._fetch(self.url)

    @DataURLMixin.safe_fetch
    def _fetch(self, url) -> pd.DataFrame:
        print("Fetching from, ", url)
        return pd.read_csv(url, sep="\t")
    
    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"PCAWG.{self.omic}.raw.parquet")
        self.read_cache = lambda cache: pd.read_parquet(cache)
        self.write_cache = lambda data, cache: data.to_parquet(cache, engine="pyarrow")

    @classmethod
    def get(cls, data: str) -> pd.DataFrame:
        return cls(PCAWG_XENA_DATASETS[data]).raw_data

    @classmethod
    def build_cache(
        cls,
        cache_dir: Path = CACHE_DIR,
        overwrite: bool = False
    ) -> None:
        for pcawg_data in PCAWG_XENA_DATASETS.values():
            print(f"Building {pcawg_data}")
            cls(pcawg_data, cache_dir).raw_data

get_pcawg = PCAWGXenaLoader.get
"""Shortcut for TCGAXenaLoader.build_cache"""

build_pcawg = PCAWGXenaLoader.build_cache
"""Shortcut for TCGAXenaLoader.build_cache"""

class PCAWGRawData:
    def __init__(self):
        for data_set in [
            "copynumber",
            "rnaseq",
            "phenotype"
        ]:
            setattr(self, data_set, get_pcawg(data_set))
    
    def __repr__(self):
        return "PCAWGRawData"

class PCAWGDataSet:
    def __init__(self, ret_union: bool = False):
        self.raw_data = PCAWGRawData()
        
        # Prepare data
        copynumber = self.configure_copynumber(self.raw_data.copynumber)
        expression = self.configure_expression(self.raw_data.rnaseq)
        phenotype = self.configure_phenotype(self.raw_data.phenotype)        
        
        if ret_union is True:
            copynumber_sample_ids = copynumber["sample"].unique()
            expression_sample_ids = expression.columns
            shared_samples = np.intersect1d(copynumber_sample_ids,
                                            expression_sample_ids,
                                            assume_unique=True)

            copynumber = copynumber.query("sample in @shared_samples")
            expression = expression.loc[:, shared_samples]
            phenotype = phenotype.query("sample in @shared_samples")
            
        self.copynumber = copynumber
        self.expression = expression
        self.phenotype = phenotype

    @staticmethod
    def configure_copynumber(data):
        return data.rename(columns={"sampleID": "sample"})

    @staticmethod
    def configure_expression(data):
        return (data
                # Parse gene IDs to remove periods
                .assign(gene_id=lambda x: \
                    x.feature.str.split(".", expand=True).get(0))
                .set_index("gene_id")
                .drop("feature", axis=1))

    @staticmethod
    def configure_phenotype(data):
        return (data
                .assign(cancer=lambda df: \
                    translate_pcawg_cancer(df.dcc_project_code))
                .explode("cancer")
                .dropna()
                .rename(columns={"icgc_specimen_id": "sample"}))

    def __repr__(self):
        return "PCAWGDataSet"

def translate_pcawg_cancer(s):
    """Harmonize project cancer types to TCGA standard.
    
    Replace with another resource, if available. Can also provide option
    to classify samples based on molecular, etc features.
    """
    PCAWG_CODES = {
        'BLCA-US': 'TCGA-BLCA', 'BRCA-US': "TCGA-BRCA", 'OV-AU': "TCGA-OV",
        'PAEN-AU': "TCGA-PDAC", 'PRAD-CA': "TCGA-PRAD", 'PRAD-US': "TCGA-PRAD",
        'RECA-EU': "TCGA-KIRC|TCGA-KIRP", 'SKCM-US': "TCGA-SKCM",
        'STAD-US': "TCGA-STAD", 'THCA-US': "TCGA-THCA", 'KIRP-US': "TCGA-KIRP",
        'LIHC-US': "TCGA-LIHC", 'PRAD-UK': "TCGA-PRAD", 'LIRI-JP': "TCGA-LIHC",
        'PBCA-DE': np.nan, 'CESC-US': "TCGA-CESC", 'PACA-AU': "TCGA-PDAC",
        'PACA-CA': "TCGA-PDAC", 'LAML-KR': "TCGA-LAML", 'COAD-US': "TCGA-COAD",
        'ESAD-UK': "TCGA-ESCA", 'LINC-JP': "TCGA-LIHC", 'LICA-FR': "TCGA-LIHC",
        'CLLE-ES': np.nan, 'HNSC-US': "TCGA-HNSC", 'EOPC-DE': "TCGA-PRAD",
        'BRCA-UK': "TCGA-BRCA", 'BOCA-UK': np.nan, 'MALY-DE': "TCGA-DLBC",
        'CMDI-UK': np.nan, 'BRCA-EU': "TCGA-BRCA", 'ORCA-IN': np.nan,
        'BTCA-SG': "TCGA-CHOL", 'SARC-US': "TCGA-SARC", 'KICH-US': "TCGA-KICH",
        'MELA-AU': "TCGA-SKCM", 'DLBC-US': "TCGA-DLBC", 'GACA-CN': "TCGA-STAD",
        'PAEN-IT': "TCGA-PDAC", 'GBM-US': "TCGA-GBM", 'KIRC-US': "TCGA-KIRC",
        'LAML-US': "TCGA-LAML", 'LGG-US': "TCGA-LGG", 'LUAD-US': "TCGA-LUAD",
        'LUSC-US': "TCGA-LUSC", 'OV-US': "TCGA-OV", 'READ-US': "TCGA-READ",
        'UCEC-US': "TCGA-UCEC",
    }
    return s.rename("pcawg_id").replace(PCAWG_CODES).str.split("|")