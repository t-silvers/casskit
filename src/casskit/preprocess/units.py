# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import FunctionTransformer

from ..descriptors import OneOf


def _abstolog21p(X):
    return np.log2(1 + X)

def _log21ptoabs(X):
    return (np.exp2(X) - 1).astype(int)

log21ptoabs = FunctionTransformer(_log21ptoabs, inverse_func=_abstolog21p)

def _abstolog2r(X):
    return np.log2(X/2)

def _log2rtoabs(X):
    return np.exp2(X) * 2

log2rtoabs = FunctionTransformer(_log2rtoabs, inverse_func=_abstolog2r)


class ToCounts(BaseEstimator, TransformerMixin):
    
    TFORM_FUNCS = {
        "log2(count+1)": log21ptoabs.fit_transform,
        "counts": lambda X: X,
        "log2(copy-number/2)": log2rtoabs.fit_transform,
    }
    
    units = OneOf(*TFORM_FUNCS.keys())
    
    def __init__(self, units: str):
        self.units = units

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return self.TFORM_FUNCS[self.units](X)