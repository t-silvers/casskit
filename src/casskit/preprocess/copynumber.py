from pathlib import Path
from typing import List
import warnings

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


__all__ = ["ThinCopynumber", "DiploidGISTIC"]


class ThinCopynumber(BaseEstimator, TransformerMixin):
    """Thin copynumber data to decrease correlations."""
    def __init__(
        self,
        drop_duplicates: bool = True,
        dups_tol_abs: int = 0,
        dups_tol_corr: float = 0.9,
        thinning: bool = True,
        thin_ar: float = 0.5,
    ) -> None:
        self.drop_duplicates = drop_duplicates
        self.dups_tol_abs = dups_tol_abs
        self.dups_tol_corr = dups_tol_corr
        self.thinning = thinning
        self.thin_ar = thin_ar
    
    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.transform(X)

    def transform(self, X, y=None):
        if self.drop_duplicates:
            X_tform = self.drop_dups(X)
        
        if self.thinning:
            X_tform = self.thin(X_tform)
        
        return X_tform

    def drop_dups(self, X):
        """Drop duplicate columns."""
        if self.dups_tol_abs:
            X = self.drop_dups_abs(X)
        
        # Drop columns with correlation greater than tolerance
        if self.dups_tol_corr:
            X = self.drop_dups_corr(X)
        
        return X
    
    def drop_dups_abs(self, X):
        """Coarsen copynumber and drop duplicate columns."""
        if self.dups_tol_abs == 0:
            return X.drop_duplicates()

        else:
            X_coarse = X.round(self.dups_tol_abs)
            return X.loc[X_coarse.drop_duplicates().index]
    
    def drop_dups_corr(self, X):
        """Drop duplicate columns with absolute correlation greater than tolerance."""
        # https://stackoverflow.com/questions/17778394/list-highest-correlation-pairs-from-a-large-correlation-matrix-in-pandas
        dups = set()
        corr_matrix = X.corr().abs()
        for i in range(len(corr_matrix.columns)):
            for j in range(i):
                if corr_matrix.iloc[i, j] > self.dups_tol_corr:
                    colname = corr_matrix.columns[i]
                    dups.add(colname)
        
        return X.drop(columns=dups)

    def thin(self, X):
        """Thin copynumber data to decrease correlations."""
        # Thin copynumber data to decrease correlations
        if self.thin_ar:
            X = self.thin_ar(X)
        
        return X


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

        elif self.units == "counts":
            X_tform = X.mask(
                lambda df: (df >= self.min_dip) & (df <= self.max_dip)
            ).fillna(2)

        return X_tform