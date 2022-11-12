import logging
from typing import Dict, List
import warnings

import pandas as pd
from pathlib import Path

from casskit.descriptors import Validator
from casskit import config

logging.basicConfig(filename=config.get_logging(), level=logging.INFO)


class PPSignal(Validator):

    def __init__(
        self,
        eval_func: str = "default",
        tol: float = 0.5,
        error_tol: float = 0.9,
        error_f: float = 0.1,
    ):
        self.eval_func = self.eval_funcs.get(eval_func)
        if error_tol < tol:
            warnings.warn("Error threshold, error_tol, is "
                          "less than warning threshold, tol.")

        self.tol = tol
        self.etol = error_tol
        self.ef = error_f

    @property
    def eval_funcs(self) -> Dict:
        return {"default": self.spearman}
    
    @staticmethod
    def spearman(df):
        return df.corr().loc["original", "transformed"]

    def _calculate_pp_signal(self, X_orig: pd.DataFrame, X_tform: pd.DataFrame):
        return (pd.concat([X_orig.assign(X_from="original"),
                           X_tform.assign(X_from="transformed")], axis=0)
                .rename_axis("sample")
                .reset_index()
                .melt(id_vars=["sample", "X_from"],
                      var_name="feature",
                      value_name="value")
                .pivot_table(index=["feature", "sample"],
                             columns="X_from",
                             values="value")
                .groupby("feature")
                .apply(lambda df: self.eval_func(df))
                .dropna())

    def validate(self, data: Dict[str, pd.DataFrame]):
        X_orig = data.get("original")
        X_tform = data.get("transformed")
        
        corr_s = self._calculate_pp_signal(X_orig, X_tform)
        
        if corr_s.empty:
            raise ValueError("Prepared data is missing feature labels.")
        
        # logging.info("Preprocessing corr - {}".format(corr_s.to_string()))
        
        # Validate individual features
        if (corr_s < self.tol).any():
            warnings.warn("Poor correlation between original and transformed data "
                          "for some features. Please check your preprocessing steps.")

        if (corr_s < self.etol).any():
            if (corr_s < self.etol).sum() > self.ef:
                raise ValueError("Correlation between original and transformed data "
                                f"is below the error threshold {self.etol} for more "
                                f" than {self.ef:.2%} of features."
                                "Please check your preprocessing steps.")

        perc_above_tol = (corr_s > self.tol).value_counts(normalize=True).loc[True]
        print(f"{perc_above_tol:.2%} of features above threshold {self.tol}.")

        # Validate averages over features
        mean_val = corr_s.mean()
        if mean_val < self.tol:
            warnings.warn("Poor correlation between original and transformed data "
                          "for some features. Please check your preprocessing steps.")

        if mean_val < self.etol:
            raise ValueError("Correlation between original and transformed data "
                            f"is below the error threshold {self.etol} for more "
                            f" than {self.ef:.2%} of features."
                            "Please check your preprocessing steps.")
            
        print(f"Processed features have average corr {mean_val:.2f} with original.")