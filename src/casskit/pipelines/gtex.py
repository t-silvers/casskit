from typing import Dict, List

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
    ToCounts,
    VariationThreshold
)


class GTEx(BaseEstimator, TransformerMixin):
    """
    https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/leafcutter/src/cluster_prepare_fastqtl.py    
    """
    
    data = PPSignal(tol=0.5, error_tol=0.9, error_f=0.2)
    
    def __init__(
        self,
        units: OneOf("log2(count+1)", "counts") = "log2(count+1)",
        lib_size: np.ndarray = None,
        genes: List[str] = None,
        min_cpm: float = 1,
        max_freq_zero: float = 0.3,
        cv2_min: float = 0.8,
    ):
        self.units = units
        self.lib_size = lib_size
        self.genes = genes
        self.min_cpm = min_cpm
        self.max_freq_zero = max_freq_zero
        self.cv2_min = cv2_min

    @property
    def gtex_preprocess(self) -> Pipeline:
        return Pipeline(
            [("As counts", ToCounts(units=self.units)),
             ("TMM", EdgeRCPM(lib_size=self.lib_size)),
             ("Protein coding", ProteinCoding(genes=self.genes)),
             ("Filter low expression", CountThreshold(self.min_cpm, self.max_freq_zero)),
             ("Filter low variance", VariationThreshold(self.cv2_min)),
             ("RINT", RINT())]
        )

    def fit(self, X, y=None):
        return self.gtex_preprocess.fit(X, y)
    
    def transform(self, X):
        self.transformed = self.gtex_preprocess.transform(X)
        self.data = {"original": X, "transformed": self.transformed}
        return self.transformed