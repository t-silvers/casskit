from dataclasses import dataclass, field
from typing import List, Optional, Union

import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from casskit.models.latent.expression_pcs import BatchModelEPCS
from casskit.pipelines.clinical import ClinicalCovariates
from casskit.pipelines.gtex import GTEx
from casskit.preprocess.mutation import AggMutations
from casskit.preprocess.units import ToCounts


@dataclass
class ModelFrame:
    variants: Union[pd.DataFrame, None]
    gene_copynumber: Union[pd.DataFrame, None]
    cnvr_copynumber: Union[pd.DataFrame, None]
    expression: Union[pd.DataFrame, None]
    phenotype: Union[pd.DataFrame, None]

@dataclass
class Phenotype:
    raw: pd.DataFrame
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        return ClinicalCovariates(
            imp_strategy="most_frequent", min_frequency=0.2
        )

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
        
        self.prepared = self.preprocessing.fit_transform(self.raw)

@dataclass
class Variants:
    raw: pd.DataFrame
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        return AggMutations(
            min_vaf=0.2, min_cases=10, binary=True
        )

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
        
        self.prepared = self.preprocessing.fit_transform(self.raw)

@dataclass
class Expression:
    raw: pd.DataFrame
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    @property
    def gtex_preprocessing(self) -> Pipeline:
        return GTEx(
            units="log2(count+1)", lib_size=None,
            genes=self.raw.columns.tolist(),
            min_cpm=1, max_freq_zero=0.3, cv2_min=0.8
        )

    def prepare(self) -> Pipeline:
        return Pipeline([
            ("GTEx preprocess", self.gtex_preprocessing),
            ("Batch effects", BatchModelEPCS("elbow")),
            ("center", StandardScaler(with_mean=True, with_std=False)),
        ])

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
        
        self.prepared = self.preprocessing.fit_transform(self.raw)

@dataclass
class GeneCopyNumber:
    raw: pd.DataFrame
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        return Pipeline([("", "passthrough")])

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
        prepared = self.preprocessing.fit_transform(self.raw)
        self.prepared = (prepared
                         .droplevel(["Chromosome", "Start", "End"])
                         .transpose())

@dataclass
class BinCopyNumber:
    raw: pd.DataFrame
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        return Pipeline([("", "passthrough")])

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
        prepared = self.preprocessing.fit_transform(self.raw)
        self.prepared = (prepared
                         .droplevel(["Chromosome", "Start", "End"])
                         .transpose())

@dataclass
class CNVRCopyNumber:
    raw: pd.DataFrame
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        return Pipeline([("", "passthrough")])

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
        prepared = self.preprocessing.fit_transform(self.raw)
        self.prepared = (prepared
                         .droplevel(["Chromosome", "Start", "End"])
                         .transpose())
