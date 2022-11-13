# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import FunctionTransformer

from ..descriptors import OneOf


def abstolog21p(X):
    return np.log2(1 + X)

@FunctionTransformer(inverse_func=lambda X: abstolog21p(X))
def log21ptoabs(X):
    return (np.exp2(X) - 1).astype(int)

def abstolog2r(X):
    return np.log2(X/2)

@FunctionTransformer(inverse_func=lambda X: abstolog2r(X))
def log2rtoabs(X):
    return np.exp2(X) * 2

class ToCounts(BaseEstimator, TransformerMixin):
    
    TFORM_FUNCS = {
        "log2(count+1)": log21ptoabs,
        "counts": lambda X: X,
        "log2(copy-number/2)": log2rtoabs,
    }
    
    units = OneOf(*TFORM_FUNCS.keys())
    
    def __init__(self, units: str):
        self.units = units

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        if self.units == "log2(count+1)":
            return self.log21ptoabs(X)
        
        elif self.units == "counts":
            return X

    @staticmethod
    def log21ptoabs(x):
        return (np.exp2(x) - 1).astype(int)
