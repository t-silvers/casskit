from typing import Dict

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.utils import validation

from .base import BatchModel


class BatchModelEPCS(BatchModel):
    """Estimate and adjust for expression PCs.
    
    Adapted from: PCAForQTL https://github.com/heatherjzhou/PCAForQTL
    
    Note:
    -----
    "we use the fully processed gene expression matrix for Colon - Transverse
    from GTEx V8 (2020) as an example. In GTEx's case, “fully processed”
    means TPM normalized, filtered, TMM normalized, and inverse normal
    transformed." - https://github.com/heatherjzhou/PCAForQTL
    
    """
    def __init__(self, select_method, ncomp: int = None, **kwargs):
        super().__init__(**kwargs)
        self.select_method = select_method
        self.ncomp = ncomp

    @property
    def ncomp_methods(self) -> Dict:
        return {
            "elbow": self.elbow,
            "buja_eyuboglu": self.buja_eyuboglu,
            "top_k": self.top_k,
        }

    def elbow(self, X):
        ncomp = self.find_elbow_custom(PCA().fit(X).explained_variance_)        
        self.pca_ = PCA(n_components=ncomp).fit(X)
        self.components_ = self.pca_.components_
        self.n_components_ = self.pca_.n_components_

    @staticmethod
    def find_elbow_custom(variance: np.ndarray) -> int:
        # TODO: Justify this method
        A = np.vstack([np.arange(len(variance)), np.ones(len(variance))]).T
        m, c = np.linalg.lstsq(A, variance, rcond=None)[0]
        resids = variance - np.dot(A, [m, c])

        i = 0
        while i in range(len(resids)):
            if resids[i] < 0:
                break
            i += 1

        return i

    @staticmethod
    def find_elbow_PCAtools(variance: np.ndarray) -> int:
        """Find the elbow in a curve.

        Note:
        ----
        Uncritical translation of R code from
        https://github.com/kevinblighe/PCAtools/blob/master/R/findElbowPoint.R

        Returns:
            int: The index of the elbow.
        """                
        # Finding distance from each point on the curve to the diagonal.
        dy = -np.diff([np.amin(variance), np.amax(variance)])
        dx = len(variance) - 1
        l2 = np.sqrt(dy ** 2 + dx ** 2)
        dy /= l2
        dx /= l2

        dy0 = variance - variance[0]
        dx0 = np.arange(len(variance))
        parallel_l2 = np.sqrt((dx0 * dx) ** 2 + (dy0 * dy) ** 2)
        normal_x = dx0 - dx * parallel_l2
        normal_y = dy0 - dy * parallel_l2
        normal_l2 = np.sqrt(normal_x ** 2 + normal_y ** 2)

        # Picking the maximum normal that lies below the line.
        # If the entire curve is above the line, we just pick the last point.
        below_line = (normal_x < 0) & (normal_y < 0)
        
        if any(below_line):
            return below_line[below_line][np.argmax(normal_l2[below_line])]
        else:
            return len(variance)

    def buja_eyuboglu(self):
        raise NotImplementedError

    def top_k(self, X) -> PCA:        
        self.pca_ = PCA(n_components=self.ncomp).fit(X)
        self.components_ = self.pca_.components_
        self.n_components_ = self.pca_.n_components_

    def filter_known_covariates():
        pass

    def fit(self, X, y=None):
        """    
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        """
        self._fit(X)
        return self

    def transform(self, X, y=None):
        return self.pca_.inverse_transform(self.pca_.transform(X))

    def _fit(self, X):
        self.ncomp_methods[self.select_method](X, **self.kwargs)

    def __repr__(self):
        return f'BatchModel(method={self.select_method})'