from pathlib import Path
from typing import List
import warnings

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


__all__ = ["AggMutations"]


class AggMutations(BaseEstimator, TransformerMixin):
    """Aggregate mutation data."""
    def __init__(
        self,
        min_vaf: float = 0.2,
        min_cases: int = 10,
        binary: bool = True,
    ):
        self.min_vaf = min_vaf
        self.min_cases = min_cases
        self.binary = binary

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.transform(X)

    def transform(self, X, y=None):
        self.feature_names = X.columns
        X_agg = self.cis_mut_agg(X)
        if (self.binary and not X_agg.empty):
            return self._binarize_mutation(X_agg)
        return X_agg

    @staticmethod
    def _binarize_mutation(df):
        return df.fillna(0).mask(lambda x: x > 0, 1).astype("category")

    @staticmethod
    def _combine_multiidx(df):
        df.columns = ['_'.join(col) for col in df.columns.values]
        return df

    def cis_mut_agg(self, X):
        """Aggregate mutation types.
        
        Note:
            Use VAF >+ 20% as cutoff for aggregation. This is a rule of thumb
            that does not have a strong theoretical basis or use ploidy and
            tumor purity.
            
            J Deng, et al. (2020) "The prognostic impact of variant allele
            frequency (VAF) in TP53 mutant patients with MDS: A systematic
            review and meta-analysis"
        
            ~Use argmax(10, .1*10) as rule of thumb for observations:~
            Use 10 as rule of thumb for observations:
                
                Long JS (1997) "Regression Models for categorical and limited dependent variables"
                RD Riley, et al. (2018) "Minimum sample size for developing a multivariable prediction model: PART II - binary and time-to-event outcomes"        
        """
        try:
            return ((X_ := self._filter_by_vaf(X))
                    .pivot_table(index=["gene", "effect"], columns="Sample_ID", values="dna_vaf")
                    .reindex(self._filter_by_cases(X_))
                    .transpose()
                    .pipe(self._combine_multiidx))
        
        except KeyError:
            warnings.warn("No passing mutations found in dataset.")
            return pd.DataFrame()

    def _filter_by_vaf(self, X):
        """Tally mutations."""
        return (X
                .query("dna_vaf >= @self.min_vaf")
                .assign(effect = lambda x: x.pop("effect").str.split(';'))
                .explode("effect"))

    def _filter_by_cases(self, X):
        """Tally mutations.
        
        More performant than comparable transformers.
        """
        return (X
                .value_counts(subset=["gene", "effect"])
                .mask(lambda x: x < self.min_cases)
                .dropna()
                .index)

