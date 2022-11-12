# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import numpy as np
import pandas as pd
import qtl
from scipy.special import erfinv
from scipy.stats import rankdata, norm
from sklearn.base import TransformerMixin
from sklearn.preprocessing import FunctionTransformer


@FunctionTransformer
def rank_inverse_normal(A, k: float = 3.0/8):
    """Rank inverse normal transform, (R)INT
    
    INT(s) = Φ^-1 [(rank(s) - k) / (n - 2k + 1)]

    where   Φ is the standard normal cumulative distribution function / probit function
            k is an adjustable offset, the Blom offset, 3/8, by default
            n is the number of observations
    
    Note:
        - https://cran.r-project.org/web/packages/RNOmni/vignettes/RNOmni.html
    """
    return norm.ppf((rankdata(A, method="average", axis=1) - k) / (A.shape[1] - 2*k + 1))


class ExpressionPreprocess(TransformerMixin):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.kwargs = kwargs
    
    @property
    def library_size(self) -> pd.Series:
        """Library size used to normalize the data."""
        if self.units == 'count':
            return self.data.sum(axis=0).astype(int)

        return self.data.sum(axis=0)





