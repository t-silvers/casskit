from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


class AggMutations(BaseEstimator, TransformerMixin):
    """Aggregate mutation data."""
    def __init__(
        self,
        min_vaf: float = 0.2,
        min_cases: int = 10
    ):
        self.min_vaf = min_vaf
        self.min_cases = min_cases

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.fit(X, y).transform(X)

    def transform(self, X, y=None):
        self.feature_names = X.columns
        return self.cis_mut_agg(X)

    def get_feature_names_out(self, input_features=None):
        return self.feature_names

    @staticmethod
    def encode_mutation(df):
        return df.mask(lambda df: df > 0, "mut").mask(lambda df: df == 0, "wt")

    def cis_mut_agg(self, X):
        """Aggregate mutation types.
        
        Note:
            Use VAF >+ 20% as cutoff for aggregation. This is a rule of thumb that does not have a strong theoretical basis
            or use ploidy and tumor purity.
                J Deng, et al. (2020) "The prognostic impact of variant allele frequency (VAF) in TP53 mutant patients with MDS: A systematic review and meta-analysis"
        
            ~Use argmax(10, .1*10) as rule of thumb for observations:~
            Use 10 as rule of thumb for observations:
                Long JS (1997) "Regression Models for categorical and limited dependent variables"
                RD Riley, et al. (2018) "Minimum sample size for developing a multivariable prediction model: PART II - binary and time-to-event outcomes"        
        """
        filtered_data = (X
                         .query("dna_vaf >= @self.min_vaf")
                         .groupby(['ensembl_gene_id','mutation_effect'])
                         .count()
                         .query("sample >= @self.min_cases"))
        
        if filtered_data.empty:
            return filtered_data
        
        else:
            second_pass = lambda df: df.loc[df.mask(lambda x: x == 'wt').count(axis=1).mask(lambda x: x < self.min_cases).index]
            return (filtered_data
                    .filter([])
                    .join(X.set_index(['ensembl_gene_id', 'mutation_effect']))
                    .reset_index()
                    .assign(mutation_eqtl = lambda x: x['ensembl_gene_id'] + '_' + x['mutation_effect'])
                    .pivot_table(index=['ensembl_gene_id', 'mutation_eqtl'], columns='sample', values='dna_vaf', fill_value=0)
                    .pipe(self.encode_mutation)
                    # Second check for minimum number of cases
                    .pipe(second_pass))
