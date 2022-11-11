# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from abc import ABC

from sklearn.base import BaseEstimator, TransformerMixin


class BatchModel(TransformerMixin, BaseEstimator, ABC):
    """Batch effects model.

    Calculate and adjust for latent variables that capture batch effects.

    Note:
        - W Mao, et al. (2020) "DataRemix: a universal data transformation for optimal inference from gene expression datasets"
        - HJ Zhou, et al. (2022) "PCA outperforms popular hidden variable inference methods for molecular QTL mapping"
        - ComBat, ComBat-Seq, sva, svaseq [sva](https://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf)

    Args:
        expression (pandas DataFrame): Expression data as sample by gene matrix.
        method (str): How to estimate latent variables. Default is 'epcs'.

    Attributes:

    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.kwargs = kwargs
