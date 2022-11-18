from dataclasses import dataclass, field
from typing import List, Optional, Union

from cneqtl.cnvrs import BinCaller, GeneCaller, CPDcaller
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from casskit import dask_cluster
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
            ("GTEx preprocess", self.gtex_preprocess),
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
    prefix: str
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        gene_caller = GeneCaller(
            assembly="GRCh37", nb_cpu=1, ddf_nparts=100, prefix=self.prefix
        )
        
        return Pipeline([
            ("gene copynumber", gene_caller),
            ("to counts", ToCounts(units="log2(copy-number/2)")),
        ])

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
    
        if not self.preprocessing["gene copynumber"].path_cache.exists():
            with dask_cluster(4, '32GB', 20, 1, "2:00:00") as (cluster, client):
                prepared = self.preprocessing.fit_transform(self.raw)
        else:
            prepared = self.preprocessing.fit_transform(self.raw)

        self.prepared = (prepared
                         .droplevel(["Chromosome", "Start", "End"])
                         .transpose())

@dataclass
class BinCopyNumber:
    raw: pd.DataFrame
    prefix: str
    bin_cnvr_width: int = int(1E5)
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        return BinCaller(
            width=self.bin_cnvr_width, ddf_nparts=200, prefix=self.prefix
        )

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
    
        if not self.preprocessing.path_cache.exists():
            with dask_cluster(4, '32GB', 20, 1, "2:00:00") as (cluster, client):
                self.prepared = self.preprocessing.fit_transform(self.raw)
        else:
            self.prepared = self.preprocessing.fit_transform(self.raw)

@dataclass
class CNVRCopyNumber:
    raw: pd.DataFrame
    prefix: str
    bin_cnvr_width: int = int(1E5)
    preprocessing: Optional[Pipeline] = field(init=True, default=None)
    prepared: pd.DataFrame = field(init=False)

    def prepare(self) -> Pipeline:
        cnvr_caller = CPDcaller(
            width=self.bin_cnvr_width, ddf_nparts=200,
            cpd_search_method="kcpd", ragged=True, prefix=self.prefix,
        )
        
        return Pipeline([
            ("call CNVRs", cnvr_caller),
            ("to counts", ToCounts(units="log2(copy-number/2)")),
        ])

    def __post_init__(self):
        if self.preprocessing is None:
            self.preprocessing = self.prepare()
    
        if not self.preprocessing["call CNVRs"].path_cache.exists():
            with dask_cluster(4, '32GB', 20, 1, "2:00:00") as (cluster, client):
                prepared = self.preprocessing.fit_transform(self.raw)
        else:
            prepared = self.preprocessing.fit_transform(self.raw)

        self.prepared = (prepared
                         .droplevel(["Chromosome", "Start", "End"])
                         .transpose())