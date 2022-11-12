# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import numpy as np
import pandas as pd
from qtl import norm as qtl_norm
from scipy.special import erfinv
from scipy.stats import rankdata, norm
from sklearn.base import TransformerMixin
from sklearn.preprocessing import BaseEstimator, FunctionTransformer

from ..descriptors import OneOf


@FunctionTransformer
def edger_cpm(A):
    return qtl_norm.edger_cpm(A.T).T

@FunctionTransformer
def rank_inverse_normal(A, k: float = 3.0/8):
    """Rank inverse normal transform, (R)INT
    
    INT(s) = Φ^-1 [(rank(s) - k) / (n - 2k + 1)]

    where   Φ is the standard normal cumulative distribution function / probit function
            k is an adjustable offset, the Blom offset, 3/8, by default
            n is the number of observations
    
    Note:
        - https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html
    """
    return norm.ppf((rankdata(A, method="average", axis=1) - k) / (A.shape[1] - 2*k + 1))

class ToCounts(BaseEstimator, TransformerMixin):
    UNITS = ["log2(count+1)", "counts"]
    units = OneOf(*UNITS)
    
    def __init__(self, units: str = "log2(count+1)"):
        self.units = units

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        if self.units == "log2(count+1)":
            return self.log21p_to_abs(X)
        elif self.units == "counts":
            return X

    @staticmethod
    def log21p_to_abs(x):
        return (np.exp2(x) - 1).astype(int)





