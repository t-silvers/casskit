from typing import List

import numpy as np
import pandas as pd
from qtl import norm as qtl_norm
from scipy import stats
from sklearn.base import BaseEstimator, TransformerMixin

from casskit.io import get_ensembl


class VariationThreshold(BaseEstimator, TransformerMixin):
    def __init__(self, cv2_min: float = 0.1):
        self.cv2_min = cv2_min

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return self.variation_threshold(X, self.cv2_min)

    @staticmethod
    def variation_threshold(A, cv2_min):
        ix = (stats.variation(A, axis=0, nan_policy="omit") >= cv2_min).nonzero()[0]
        return np.take(A, ix, axis=1)

class EdgeRCPM(BaseEstimator, TransformerMixin):
    def __init__(self, lib_size: np.ndarray = None):
        self.lib_size = lib_size

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return self.edger_cpm(X.T, self.lib_size).T

    @staticmethod
    def edger_cpm(counts_df, lib_size=None):
        """
        Reproduces qtl_norm.edger_cpm, which reproduces edgeR::cpm.DGEList,
        with argument lib_size.
        """
        if lib_size is None:
            lib_size = counts_df.sum(axis=0)
        tmm = qtl_norm.edger_calcnormfactors(counts_df)
        lib_size *= tmm
        
        return counts_df / lib_size * 1e6

class RINT(BaseEstimator, TransformerMixin):
    def __init__(self, k: float = 3.0/8):
        self.k = k

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return pd.DataFrame(self.rank_inverse_normal_transform(X, self.k),
                            index=X.index, columns=X.columns)

    @staticmethod
    def rank_inverse_normal_transform(A, k: float = 3.0/8):
        """Rank inverse normal transform, (R)INT
        
        INT(s) = Φ^-1 [(rank(s) - k) / (n - 2k + 1)]

        where   Φ is the standard normal cumulative distribution function / probit function
                k is an adjustable offset, the Blom offset, 3/8, by default
                n is the number of observations
        
        Note:
            - https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html
        """
        return stats.norm.ppf((stats.rankdata(A, method="average", axis=0) - k) / (A.shape[0] - 2*k + 1))

class CountThreshold(BaseEstimator, TransformerMixin):
    def __init__(self, min_cpm: int = 1, max_freq_zero: float = 0.3):
        self.min_cpm = min_cpm
        self.max_freq_zero = max_freq_zero

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return self.filter_by_count(X, self.min_cpm, self.max_freq_zero)

    @staticmethod
    def filter_by_count(A, k, n):
        """
        This function implements the filtering strategy that was intuitively
        described by Chen et al (2016). Roughly speaking, the strategy keeps genes
        that have at least min.count reads in a worthwhile number samples. More
        precisely, the filtering keeps genes that have count-per-million (CPM) above
        k in n samples, where k is determined by min.count and by the sample library
        sizes and n is determined by the design matrix.

        n is essentially the smallest group sample size or, more generally, the
        minimum inverse leverage of any fitted value. If all the group sizes are
        larger than large.n, then this is relaxed slightly, but with n always
        greater than min.prop of the smallest group size (70% by default).

        In addition, each kept gene is required to have at least min.total.count
        reads across all the samples.
        """

        ix = ((np.count_nonzero(A >= k, axis=0) / A.shape[0]) >= n).nonzero()[0]
        return np.take(A, ix, axis=1)

class ProteinCoding(BaseEstimator, TransformerMixin):
    """Filter out non-protein coding genes
    
    If in pipelin, must be run before other filtering steps.
    
    Parameters
    ----------
    genes: List[str]
        List of gene names (as Ensembl gene ID) to filter
    """
    def __init__(self, genes, assembly: str = "GRCh37"):
        self.genes = genes
        self.assembly = assembly
        self.protein_coding_genes = self.get_protein_coding(assembly)

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        ix = self.filter_protein_coding(self.genes, self.protein_coding_genes)
        return np.take(X, ix, axis=1)

    @staticmethod
    def get_protein_coding(assembly) -> List:
        return (get_ensembl(assembly)
                .query("gene_biotype == 'protein_coding'")
                ["gene_id"]
                .tolist())
    
    @staticmethod
    def filter_protein_coding(X_genes, protein_coding_genes):
        protein_coding = set(X_genes).intersection(protein_coding_genes)
        return np.sort([X_genes.index(x) for x in protein_coding])