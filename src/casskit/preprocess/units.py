import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import FunctionTransformer

from casskit.descriptors import OneOf


__all__ = ["UNITS", "ToCounts"]


UNITS = {
    "log1p": "log1p",
    "log2": "log2",
    "log10": "log10",
    "log2(count+1)": "log1p",
    "log2(copy-number/2)": "",
    "absolute": "counts",
    "abs": "counts",
    "counts": "counts"
}

# TODO: FunctionTransformer arg check_inverse doesn't work when input has no dtype (eg pd.DataFrame)


def _abstolog21p(X):
    return np.log2(1 + X)

def _log21ptoabs(X):
    return (np.exp2(X) - 1).astype(int)

log21ptoabs = FunctionTransformer(
    _log21ptoabs, inverse_func=_abstolog21p, check_inverse=False
)

def _abstolog2r(X):
    return np.log2(X/2)

def _log2rtoabs(X):
    return np.exp2(X) * 2

log2rtoabs = FunctionTransformer(
    _log2rtoabs, inverse_func=_abstolog2r, check_inverse=False
)


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