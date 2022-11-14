from typing import Dict, List, Tuple, Union

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer, make_column_selector
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

from casskit.preprocess.phenotype import ReplaceNA, TumorStageEncoder


class ClinicalCovariates(BaseEstimator, TransformerMixin):

    TCGA_VARS = [
        "batch_number", "ethnicity.demographic", "gender.demographic",
        "race.demographic", "age_at_diagnosis.diagnoses", "clinical_stage",
        "tumor_stage.diagnoses", "neoplasm_histologic_grade"
    ]
    
    CAT_VARS = [
        "batch_number", "ethnicity.demographic",
        "gender.demographic", "race.demographic",
    ]

    STAGE_VARS = ["clinical_stage", "tumor_stage.diagnoses"]
    
    def __init__(
        self,
        imp_strategy: str = "most_frequent",
        onehot_drop: str = "if_binary",
        min_frequency: Union[int, float] = 0.2,
    ):
        self.imp_strategy = imp_strategy
        self.onehot_drop = onehot_drop
        self.min_frequency = min_frequency

    def fit(self, X, y=None):
        self.clinical_preprocessor.fit(X[self.TCGA_VARS])
        self.fitted_ = True
        
        return self

    def transform(self, X):
        if not self.fitted_:
            self.fit(X[self.TCGA_VARS])

        return self.clinical_preprocessor.transform(X[self.TCGA_VARS])
    
    def fit_transform(self, X, y=None):
        """Custom fit_transform method for checks."""
        return self.clinical_preprocessor.fit_transform(X[self.TCGA_VARS])
    
    @property
    def clinical_preprocessor(self) -> Pipeline:
        return Pipeline(
            [("null_vals", self.null_values),
             ("stage_encoder", self.stage_encoder),
             ("impute1", self.impute),
             ("one_hot", self.one_hot),
             ("var_thresh", self.singular)]
        )
    
    @property
    def null_values(self) -> ColumnTransformer:
        return ColumnTransformer(transformers=[
            ("null", ReplaceNA(),
             make_column_selector(dtype_include=[object, np.number]))
            ], verbose_feature_names_out=False, remainder="passthrough")

    @property
    def impute(self) -> ColumnTransformer:
        return ColumnTransformer(transformers=[
            ("impute", SimpleImputer(strategy=self.imp_strategy),
             make_column_selector(dtype_include=[object, np.number]))
            ], verbose_feature_names_out=False, remainder="passthrough")

    @property
    def one_hot(self) -> ColumnTransformer:
        return ColumnTransformer(transformers=[
            ("onehot", OneHotEncoder(
                min_frequency=self.min_frequency,
                drop=self.onehot_drop,
                sparse_output=False    
            ), self.CAT_VARS)
            ], verbose_feature_names_out=False, remainder="passthrough")

    @property
    def stage_encoder(self) -> ColumnTransformer:
        return ColumnTransformer(transformers=[
            ("stage", TumorStageEncoder(), self.STAGE_VARS)
            ], verbose_feature_names_out=False, remainder="passthrough")

    @property
    def singular(self) -> ColumnTransformer:
        return ColumnTransformer(transformers=[
            ("singular", VarianceThreshold(threshold=0),
             make_column_selector(dtype_include=np.number))
            ], verbose_feature_names_out=False, remainder="passthrough")
