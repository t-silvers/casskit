from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


class ReplaceNA(BaseEstimator, TransformerMixin):
    """Replace non-standard missing values with np.nan."""
    def __init__(self, replace_vals: List[str] = ["not reported"]):
        self.replace_vals = replace_vals

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.fit(X, y).transform(X)

    def transform(self, X, y=None):
        self.feature_names = X.columns
        X = X.copy()

        for rval in self.replace_vals:
            mask = X.apply(lambda x: x.astype(str).str.lower()) == rval
            X[mask] = np.nan

        return X

    def get_feature_names_out(self, input_features=None):
        return self.feature_names




class CisMutModel(SampleModel):
    
    def __init__(self, path: Path, expression: pd.DataFrame, mutation: pd.DataFrame, covariates: pd.DataFrame = None, **kwargs):
        self.mutation = mutation
        self.mut_method = kwargs.get('mut_method', 'agg')
        
        super().__init__(path=path, expression=expression, covariates=covariates, **kwargs)
        self.mutation_cache = Path(self.mut_directory, 'mutation_eqtls.feather')

    def _remove_mut_model(self):
        self.mutation_cache.unlink(missing_ok=True)

    @property
    def mutation_eqtls(self):
        return self.fit_mut_model()

    @staticmethod
    def cis_mut_agg(data: pd.DataFrame, min_vaf=0.2, min_cases=10):
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
        # Minimum number of cases for each mutation type
        # min_cases = max(min_cases, int(.1 * data.shape[0]))
        min_cases = 10
        
        encode_mutation = lambda df: df.mask(lambda df: df > 0, "mut").mask(lambda df: df == 0, "wt")
        filtered_data = (data.query("dna_vaf >= @min_vaf")
                         .groupby(['ensembl_gene_id','mutation_effect'])
                         .count()
                         .query("sample >= @min_cases"))
        
        if filtered_data.empty:
            return filtered_data
        
        else:
            second_pass = lambda df: df.loc[df.mask(lambda x: x == 'wt').count(axis=1).mask(lambda x: x < min_cases).index]
            return (filtered_data
                    .filter([])
                    .join(data.set_index(['ensembl_gene_id', 'mutation_effect']))
                    .reset_index()
                    .assign(mutation_eqtl = lambda x: x['ensembl_gene_id'] + '_' + x['mutation_effect'])
                    .pivot_table(index=['ensembl_gene_id', 'mutation_eqtl'], columns='sample', values='dna_vaf', fill_value=0)
                    .pipe(encode_mutation)
                    # Second check for minimum number of cases
                    .pipe(second_pass))

    def fit_mut_model(self) -> pd.DataFrame:
        """Mutation eQTLs model.

        Note:
            aggregation eQTL (aeQTL) is a method to identify cis-eQTLs from mutation data.

        Args:
            expression (pandas DataFrame): Expression data as sample by gene matrix.
            method (str): How to map mutation eQTLs. Default is 'agg'.

        Returns:
            pandas DataFrame: Sample by [ensembl_gene_id, mutation_eqtl] eQTL matrix.
        """
        if not self.mutation_cache.exists():
            print("Mapping mutation eQTLs...")

            if self.mut_method == 'aeqtl':
                expression_data = self.sample_residualized_expression
                raise ValueError('mut_method not supported yet')
            
            elif self.mut_method == 'agg':
                mutation_eqtls = self.cis_mut_agg(self.mutation)

            mutation_eqtls.reset_index().to_feather(self.mutation_cache)

        if pd.read_feather(self.mutation_cache).empty:
            return pd.read_feather(self.mutation_cache)
        return pd.read_feather(self.mutation_cache).set_index(['ensembl_gene_id', 'mutation_eqtl']).transpose()
