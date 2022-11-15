from sklearn.base import clone
from sklearn.compose import TransformedTargetRegressor

__all__ = ["TransformedLHSRegressor"]


class TransformedLHSRegressor(TransformedTargetRegressor):
    """Wrapper for TransformedTargetRegressor to allow y as non-np.ndarray."""
    def __init__(self, regressor, transformer):
        super().__init__(
            regressor, transformer=transformer, check_inverse=False
        )

    def fit(self, X, y, **fit_params):
        self._fit_transformer(y)
        y_tform = self.transformer_.transform(y)
        
        self.regressor_ = clone(self.regressor)
        self.regressor_.fit(X, y_tform, **fit_params)

        if hasattr(self.regressor_, "feature_names_in_"):
            self.feature_names_in_ = self.regressor_.feature_names_in_

        return self