from ast import Return
from dataclasses import dataclass, field
from typing import List, Optional, Union
import warnings

import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from casskit.models.latent.expression_pcs import BatchModelEPCS
from casskit.pipelines.clinical import ClinicalCovariates
from casskit.pipelines.gtex import GTEx
from casskit.preprocess.mutation import AggMutations
from casskit.preprocess.units import ToCounts
from casskit.typing import DATAFRAME


@dataclass(frozen=False)
class ModelFrame:
    model_frame: DATAFRAME = field(default=None)

    index: List[str] = None
    impute: Union[bool, str] = True
    impute_method: str = "most_frequent"
    
    def load(self, **kwargs):
        self.expression = kwargs.get("expression", None)

        self.cnvr_copynumber = kwargs.get("cnvr_copynumber", None)
        self.gene_copynumber = kwargs.get("gene_copynumber", None)
        self.phenotype = kwargs.get("phenotype", None)
        self.variants = kwargs.get("variants", None)

        self.index = kwargs.get("index")
        self.impute = kwargs.get("impute", True)
        self.impute_method = kwargs.get("impute_method", "most_frequent")
        
        self.model_frame = self._make_model_frame()
        
        return self
    
    def _make_model_frame(self):
        # TODO: Add checks for column names
        # TODO: Add prefixes (based on arg?)

        model_frame = pd.concat([
            self.expression.add_prefix("expr_"),
            self.cnvr_copynumber,
            self.gene_copynumber.add_prefix("cn_"),
            self.phenotype,
            self.variants,
        ], axis=1)
        
        # Indices of components
        self.indices = {}
        extent = 0
        for k in [
            "expression", "cnvr_copynumber", "gene_copynumber", "phenotype", "variants"
        ]:
            nvars = getattr(self, k).shape[1]
            self.indices[k] = range(extent, nvars)
            extent += nvars
        
        model_frame = self._process_model_frame(model_frame)
        self._validate_model_frame(model_frame)
        
        return model_frame

    def _process_model_frame(self, model_frame):
        if self.index is not None:
            model_frame = model_frame.reindex(self.index)
        
        if self.impute:
            if self.impute_method == "most_frequent":
                imp = SimpleImputer(strategy="most_frequent")
            elif self.impute_method == "mean":
                imp = SimpleImputer(strategy="mean")
            else:
                raise ValueError(f"impute_method must be 'most_frequent' ",
                                 f"or 'mean' not {self.impute_method}")

            model_frame = imp.fit_transform(model_frame)

        else:
            if self.impute_method == "override":
                pass
            else:
                warnings.warn("Impute is False, so rows with NaNs are dropped. "
                            "To keep rows with NaNs, set impute='override'.")
                model_frame = model_frame.dropna()
    
        return model_frame

    def _validate_model_frame(self, model_frame):
        self.N, self.P = model_frame.shape

    def __post_init__(self):
        if self.model_frame is not None:
            self.model_frame = self._process_model_frame(self.model_frame)
            self._validate_model_frame(self.model_frame)

    def __getitem__(self, key):
        return self.model_frame[key]

    def __repr__(self):
        if self.model_frame is not None:
            return f"ModelFrame with {self.N} samples and {self.P} variables"
        return "ModelFrame with no data"

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
