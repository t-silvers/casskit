from distutils.log import warn
from logging import warning
from typing import Dict, List
import warnings

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline

from casskit.descriptors import OneOf
from casskit.preprocess.utils import PPSignal
from casskit.preprocess.expression import (
    CountThreshold,
    EdgeRCPM,
    ProteinCoding,
    RINT,
    VariationThreshold
)
from casskit.preprocess.units import ToCounts


class GTEx(BaseEstimator, TransformerMixin):
    """
    https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/leafcutter/src/cluster_prepare_fastqtl.py    
    """
    
    def __init__(
        self,
        units: OneOf("log2(count+1)", "counts") = "log2(count+1)",
        lib_size: np.ndarray = None,
        rint: bool = True,
        genes: List[str] = None,
        min_cpm: float = 1,
        max_freq_zero: float = 0.3,
        cv2_min: float = 0.8,
        validate: bool = True,
        pp_tol: float = 0.9,
        pp_etol: float = 0.8,
        pp_efrac: float = 0.2,
    ):
        self.units = units
        self.lib_size = lib_size
        self.rint = rint
        self.genes = genes
        self.min_cpm = min_cpm
        self.max_freq_zero = max_freq_zero
        self.cv2_min = cv2_min
        self.validate = validate
        self.pp_tol = pp_tol
        self.pp_etol = pp_etol
        self.pp_efrac = pp_efrac

    @property
    def validator(self) -> None:
        return PPSignal(tol=self.pp_tol, error_tol=self.pp_etol, error_f=self.pp_efrac)

    @property
    def gtex_preprocess(self) -> Pipeline:
        gtex_pipe = [
            ("As counts", ToCounts(units=self.units)),
            ("TMM", EdgeRCPM(lib_size=self.lib_size)),
            ("Protein coding", ProteinCoding(genes=self.genes)),
            ("Filter low expression", CountThreshold(self.min_cpm, self.max_freq_zero)),
            ("Filter low variance", VariationThreshold(self.cv2_min))
        ]
        
        if self.rint is True:
            gtex_pipe.append(("RINT", RINT()))

        return Pipeline(gtex_pipe)

    def fit(self, X, y=None):
        self.gtex_preprocess.fit(X)
        return self

    def transform(self, X):
        self.fit(X)
        return self.gtex_preprocess.transform(X)
    
    def fit_transform(self, X, y=None):
        """Custom fit_transform method for checks."""
        self.transformed = self.gtex_preprocess.fit_transform(X)
        try:
            self.validator.validate({
                "original": X, "transformed": self.transformed, "type": "expression"
            })
        except ValueError as e:
            if self.validate is True:
                raise e
            else:
                # Ignore validation errors
                warnings.warn(e)
                pass

        return self.transformed