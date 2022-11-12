# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from typing import Dict

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline

from ..descriptors import OneOf
from ..preprocess import EdgeRCPM, RINT, ToCounts


class GTEx(BaseEstimator, TransformerMixin):
    """
    https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/leafcutter/src/cluster_prepare_fastqtl.py
    
    we use the fully processed gene expression matrix for Colon - Transverse
    from GTEx V8 (2020) as an example. In GTEx's case, “fully processed”
    means TPM normalized, filtered, TMM normalized, and inverse normal
    transformed
    """
    def __init__(self, units: OneOf("log2(count+1)", "counts") = "log2(count+1)"):
        self.units = units

    @property
    def gtex_preprocess(self) -> Pipeline:
        return Pipeline([
            ("to counts", ToCounts(units=self.units)),
            ("TMM", EdgeRCPM),
            ("RINT", RINT)
        ])

    def fit(self, X, y=None):
        return self.gtex_preprocess.fit(X, y)
    
    def transform(self, X, y=None):
        return self.gtex_preprocess.fit(X, y)