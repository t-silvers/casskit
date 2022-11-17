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
        frac_test: float = .1,
        corr_window: int = 100,
        tol: float = 0.8,
        error_tol: float = 0.2,
        error_f: float = 0.1,
    ):
        self.eval_func = self.eval_funcs.get(eval_func)
        if error_tol > tol:
            warnings.warn("Error threshold, error_tol, is "
                          "greater than warning threshold, tol.")

        self.frac_test = frac_test
        self.corr_window = corr_window
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

    def _cnvrbins_pp_signal(self, X_orig: pd.DataFrame, X_tform: pd.DataFrame):
        return (X_tform
                .sample(frac=self.frac_test, replace=False)
                .melt(ignore_index=False)
                .reset_index()
                .merge(X_orig, suffixes=('_bin', '_seg'), on=['sample', 'Chromosome'], how="left")
                .query("(Start_bin >= Start_seg) & (End_bin <= End_seg)")
                .sort_values(["Chromosome", "Start_seg", "End_seg"])
                .filter(like="value_")
                .rolling(self.corr_window)
                .corr(pairwise=True)
                .rename_axis(["ix", "var"])
                .query("var == 'value_bin'")
                .loc[:, 'value_seg']
                .dropna()
                .reset_index(drop=True))

    def _cnvrmvcpd_pp_signal(self, X_orig: pd.DataFrame, X_tform: pd.DataFrame):
        X_orig = (X_orig.melt(var_name="sample", ignore_index=False).reset_index().astype({"Start": int, "End": int}))

        return (X_tform
                .sample(frac=self.frac_test, replace=False)
                .melt(var_name="sample", ignore_index=False)
                .reset_index()
                .astype({"Start": int, "End": int})
                .merge(X_orig, suffixes=("_tformed", "_original"), on=['sample', 'Chromosome'], how="left")
                .query("(Start_original >= Start_tformed) & (End_original <= End_tformed)")
                .groupby("Chromosome")
                .apply(
                    lambda df: (
                        df
                        .set_index("cnvr_id")
                        .sort_values("Start_tformed")
                        .filter(like="value_")
                        .rolling(100)
                        .corr(pairwise=True)
                        .dropna()
                        .rename_axis(["cnvr_id", "var"])
                        .query("var == 'value_original'")
                        .loc[:, 'value_tformed']
                        .reset_index(drop=True)
                    )))

    def validate(self, data: Dict[str, pd.DataFrame]):
        X_orig = data.get("original")
        X_tform = data.get("transformed")
        data_type = data.get("type")
        
        # TODO: Refactor to abstract
        if data_type in ["copynumber_gene", "copynumber_bins"]:
            corr_s = self._cnvrbins_pp_signal(X_orig, X_tform)

        elif data_type == "copynumber_mvcpd":
            corr_s = self._cnvrmvcpd_pp_signal(X_orig, X_tform)
        
        elif data_type == "expression":
            corr_s = self._calculate_pp_signal(X_orig, X_tform)
        
        if corr_s.empty:
            raise ValueError("Prepared data is missing feature labels.")
        
        # logging.info("Preprocessing corr - {}".format(corr_s.to_string()))
        
        # Validate individual features
        if (corr_s < self.tol).any():
            warnings.warn("Poor correlation between original and transformed data "
                          "for some features. Please check your preprocessing steps.")

        if (corr_s < self.etol).any():
            try:
                if (corr_s < .8).value_counts(normalize=True).loc[True] > self.ef:
                    raise ValueError("Correlation between original and transformed data "
                                    f"is below the error threshold {self.etol} for more "
                                    f" than {self.ef:.2%} of features."
                                    "Please check your preprocessing steps.")
            except KeyError:
                raise KeyError("Correlation between original and transformed data "
                               "is below the error threshold for all features.")

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
                            f"than {self.ef:.2%} of features."
                            "Please check your preprocessing steps.")
            
        print(f"Processed features have average corr {mean_val:.2f} with original.")