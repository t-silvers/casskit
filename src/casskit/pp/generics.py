import numpy as np
from scipy import stats
from sklearn.base import BaseEstimator, TransformerMixin


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
