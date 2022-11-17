# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from pathlib import Path
from typing import Optional

import pandas as pd
from scipy.stats import norm
from sklearn.experimental import enable_iterative_imputer # explicitly require experimental feature
from sklearn.impute import IterativeImputer, KNNImputer

import casskit.io.base as base
import casskit.utils as utils
import casskit.config as config


class TCGATumorPurityAran2015(base.ElsevierLink):
    """Tumor purity estimates.
    
    Notes
    -----
    Systematic pan-cancer analysis of tumour purity. Aran et al. 2015.
    DOI: https://doi.org/10.1038/ncomms9971
    
    #TODO: Why are there values in column 7?
    CPE = Consensus Purity Estimate

    """
    
    URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"
    
    def __init__(
        self,
        cache_dir: Optional[Path] = None,
        impute: bool = True,
        impute_method: str = "knn", # "iterative" broken in dev branch
        ret_recommended: bool = True,
    ):
        if cache_dir is None:
            cache_dir = Path(config.CACHE_DIR)

        super().__init__(
            url=self.URL,
            skiprows=3,
            cache_name="tumor_purity_aran2015",
            cache_dir=cache_dir,
        )

        self.impute = impute # Not implemented
        self.impute_method = impute_method # Not implemented
        self.ret_recommended = ret_recommended        
        self.omic = "Consensus Purity Estimate"
        self.units = "proportion"

    @property
    def prepared_data(self) -> pd.DataFrame:
        return self.prepare()

    def impute_missing(self, data: pd.DataFrame) -> pd.DataFrame:
        if self.impute_method == "iterative":
            # TODO: IterativeImputer is broken in dev branch with output="pandas"
            # TODO: Performance https://scikit-learn.org/stable/auto_examples/impute/plot_missing_values.html#missing-information
            imp = IterativeImputer(imputation_order="roman")

            # Order of columns is important for IterativeImputer?
            col_order = data.count().drop("cpe").sort_values(ascending=False).index.tolist() + ["cpe"]

            data_imputed = (data.copy()
                            # Probit transform. Subtract constant to avoid inf.
                            .sub(1E-3)
                            .apply(norm.ppf)
                            .filter(col_order)
                            .pipe(imp.fit_transform))
            
            # Inverse probit transform and return as df
            return pd.DataFrame(norm.cdf(data_imputed), index=data.index, columns=data.columns)

        elif self.impute_method == "knn":
            imp = KNNImputer(n_neighbors=10)
            return imp.fit_transform(data)

    def return_recommendation(self, data: pd.DataFrame) -> pd.DataFrame:
        """Return the CPE column."""
        return data["cpe"].rename("tumor_purity").to_frame() if self.ret_recommended is True else data
        
    def prepare(self) -> pd.DataFrame:
        return (self.raw_data
                # .drop("Unnamed: 7", axis=1)
                .iloc[:, 0:7]
                # TODO: Remove clean columns from omics class and use here
                .rename(columns={
                    "Sample ID": "sample",
                    "Cancer type": "cancer"
                })
                .pipe(utils.column_janitor)
                .set_index(["cancer", "sample"])
                .pipe(self.impute_missing)
                .pipe(self.return_recommendation))

    @classmethod
    def get_data(cls, **kwargs) -> pd.DataFrame:
        return cls(**kwargs).prepared_data

get_tumor_purity = TCGATumorPurityAran2015.get_data
"""Convenience function for tumor purity estimates from Aran et al. 2015."""
