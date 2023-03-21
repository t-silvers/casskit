from typing import List
import warnings

import numpy as np
import pandas as pd
from qtl import norm as qtl_norm
from scipy import stats
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline

from .generics import VariationThreshold
from .units import ToCounts
from ..data.io.annot import get_ensembl
from ..data.io.descriptors import OneOf


class GTEx(BaseEstimator, TransformerMixin):
    """
    https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/leafcutter/src/cluster_prepare_fastqtl.py    
    """
    
    def __init__(
        self,
        units: OneOf("log2(count+1)", "counts") = "log2(count+1)",
        lib_size: np.ndarray = None,
        rint: bool = True,
        min_cpm: float = 1,
        max_freq_zero: float = 0.3,
        cv2_min: float = 0.8,
    ):
        self.units = units
        self.lib_size = lib_size
        self.rint = rint
        self.min_cpm = min_cpm
        self.max_freq_zero = max_freq_zero
        self.cv2_min = cv2_min

    @property
    def gtex_preprocess(self) -> Pipeline:
        gtex_pipe = [
            ("As counts", ToCounts(units=self.units)),
            ("TMM", EdgeRCPM(lib_size=self.lib_size)),
            ("Protein coding", ProteinCoding()),
            ("Filter low expression", CountThreshold(self.min_cpm, self.max_freq_zero)),
            ("Filter low variance", VariationThreshold(self.cv2_min))
        ]
        
        if self.rint is True:
            gtex_pipe.append(("RINT", RINT()))

        return Pipeline(gtex_pipe)

    def fit(self, X, y=None):
        self.gtex_preprocess.fit(X)
        return self

    def transform(self, X):
        self.fit(X)
        return self.gtex_preprocess.transform(X)
    
    def fit_transform(self, X, y=None):
        """Custom fit_transform method for checks."""
        self.transformed = self.gtex_preprocess.fit_transform(X)
        return self.transformed

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
        n_obs = A.shape[0]
        a_ranked = stats.rankdata(A, method="average", axis=0)
        a_ranked_ = (a_ranked - k) / (n_obs - 2*k + 1)
        a_rint = stats.norm.ppf(a_ranked_)
        
        return a_rint

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
        freq_ge_k = (np.count_nonzero(A >= k, axis=0) / A.shape[0])
        filt_ix = (freq_ge_k >= n).nonzero()[0]
        filt_a = np.take(A, filt_ix, axis=1)
        
        return filt_a

class ProteinCoding(BaseEstimator, TransformerMixin):
    """Filter out non-protein coding genes
    
    If in pipelin, must be run before other filtering steps.
    
    Parameters
    ----------
    assembly: str
    """
    def __init__(self, assembly: str = "GRCh37"):
        self.assembly = assembly
        self.protein_coding_genes = (
            get_ensembl(assembly)
            .query("gene_biotype == 'protein_coding'")
            ["gene_id"]
            .tolist()
        )

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X_genes = X.columns
        protein_coding = list(set(X_genes).intersection(self.protein_coding_genes))
        X_protein_coding = X.reindex(columns=protein_coding)
        return X_protein_coding

    def fit_transform(self, X, y=None):
        """Custom fit_transform method for checks."""
        return self.transform(X)
