from pathlib import Path
from typing import List
import warnings

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


class ThinCopynumber(BaseEstimator, TransformerMixin):
    """Thin copynumber data to decrease correlation."""
    pass

class DiploidGISTIC(BaseEstimator, TransformerMixin):
    """Bin near-diploid values to diploid."""
    def __init__(
        self,
        min_dip: float = 1.7,
        max_dip: float = 2.3,
        units: str = "log2(copy-number/2)",
    ):
        self.min_dip = min_dip
        self.max_dip = max_dip
        self.units = units

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.transform(X)

    def transform(self, X, y=None):
        if self.units == "log2(copy-number/2)":
            X_tform = X.mask(
                lambda df: (df >= np.log2(self.min_dip/2)) & (df <= np.log2(self.max_dip/2))
            ).fillna(0)

        elif self.units == "absolute":
            X_tform = X.mask(
                lambda df: (df >= self.min_dip) & (df <= self.max_dip)
            ).fillna(2)

        return X_tform