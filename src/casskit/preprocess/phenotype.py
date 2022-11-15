from typing import List

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import OrdinalEncoder


class SimpleNaNImputer(BaseEstimator, TransformerMixin):
    """Replace non-standard missing values with np.nan.
    
    Similar behavior to
    
    >>> SimpleImputer(missing_values=replace_val, strategy="constant", fill_value=np.nan)
    
    except:
    1. Can accept a *list* of values to replace
    2. Doesn't type check input or output (i.e. can be used with data containing NaNs),
        (see sklearn.utils.assert_all_finite)
    
    """
    # can use SimpleImputer(missing_values=replace_vals, strategy="constant", fill_value=np.nan)
    def __init__(
        self, replace_vals: List[str] = ["not reported"]
    ):
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

class TumorStageEncoder(TransformerMixin, BaseEstimator):
    """Clean tumor stage data."""
    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.fit(X, y).transform(X)

    def transform(self, X):
        self.feature_names = X.columns
        
        if X.shape[1] > 1:
            X = self._combine_stage_vars(X)

        cleaned = (X.squeeze()
                   .str.upper()
                   .str.extract("STAGE ([IVX]+)", expand=False)
                   .apply(self.rom_to_int)
                   .to_frame()
                   .values)
        
        return OrdinalEncoder(
            handle_unknown="use_encoded_value", unknown_value=np.nan
        ).fit_transform(cleaned)

    @staticmethod
    def _combine_stage_vars(X):
        """If two columns, pick one with most info."""
        selected_col = (X.apply(
            lambda x: x.str.lower().str.contains("stage")
        ).sum().astype(int).idxmax())
        
        return X[selected_col]
   
    @staticmethod
    def rom_to_int(a_rom):
        """Convert a stage roman numeral to an integer.

        Adapted from https://www.oreilly.com/library/view/python-cookbook/0596001673/ch03s24.html
        Module: Roman Numerals
        Credit: Paul M. Winkler
        
        Does not handle edge cases like "I/II NOS".
        """
        if isinstance(a_rom, type(None)):
            return np.nan
        
        if not isinstance(a_rom, type("")):
            raise TypeError(f"expected string, got {a_rom}")
        
        nums = {"X":10, "V":5, "I":1}
        a_int = 0
        for i in range(len(a_rom)):
            try:
                value = nums[a_rom[i]]
                if i+1 < len(a_rom) and nums[a_rom[i+1]] > value:
                    a_int -= value
                else:
                    a_int += value

            except KeyError:
                raise ValueError(f"Input is not a valid Roman numeral: {a_rom}")
        return a_int
    
    def get_feature_names_out(self, input_features=None, use_name=False):
        if use_name:
            return np.array(self.feature_names)
        return np.array(["tumor_stage"])
