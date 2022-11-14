from typing import List

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


class TumorStageEncoder(TransformerMixin, BaseEstimator):
    def __init__(self):
        pass

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        return super().fit_transform(X, y)

    def transform(self, X):
        return self.clean_stage(X, self.cv2_min)

    def clean_stage():
        pass





import numpy as np
import pandas as pd
from pathlib import Path
from qtl import norm as qtl_norm
import scipy
from typing import Dict, List
import warnings

from ..src import utils
from ..data.tcga import TCGAdataXenaExpression, TCGAfromCache
from ..data.misc import TCGAancestryPCs


# class MappingData:
#     def __init__(self, cache_dir: Path = Path("resources"), **kwargs) -> None:
#         self.cache_directory = Path(cache_dir) / "mapping_data"
#         if not self.cache_directory.exists():
#             self.cache_directory.mkdir(parents=True)


class TCGAdataXenaPhenotype(TCGAdataXena):
    def __init__(self, cancer, primary_only=True, **kwargs):
        self.cancer = cancer
        self.omic = 'GDC_phenotype'
        self.compression='gzip'
        super().__init__(cancer, omic=self.omic, compression=self.compression, **kwargs)

        # Prepare cache if it doesn't exist
        if not self.prepared_cache.exists():
            self.raw_data = self.read_cache()
            data = self.prepare_tcga_phenotype_data(self.raw_data)
            self.save_func(self._validate_tcga_phenotype_data(data), self.prepared_cache)

    @property
    def prepared_tcga_phenotype_data(self) -> pd.DataFrame:
        return super().__call__()

    @staticmethod
    def _validate_tcga_phenotype_data(data) -> pd.DataFrame:
        return data

    @staticmethod
    def filter_cols(dat):
        """Filter out columns that are not useful for analysis."""
        INCLUDED_VARS = [
            'submitter_id.samples',
            'sample_type.samples',
            'batch_number',
            'ethnicity.demographic',
            'gender.demographic',
            'race.demographic',
            'age_at_diagnosis.diagnoses',
            'clinical_stage',
            'tumor_stage.diagnoses',
            'neoplasm_histologic_grade'
        ]
        
        return dat.filter(INCLUDED_VARS)

    @staticmethod
    def _impute_stage(data):
        """Impute stage based on clinical_stage and tumor_stage.
        Not run.
        """
        imp_mean = IterativeImputer(random_state=0)
        return imp_mean.fit_transform(data) # dat must be numeric, see below

    @staticmethod
    def _parse_stage_col(stage_col: pd.Series) -> List:
        """Parse stage column."""
        # TODO: Should just count the number of 'i' in string + something to handle errors better
        STAGE_DICT = {
            'stage i': 1, 'stage ia': 1, 'stage ia1': 1, 'stage ia2': 1, 'stage ib': 1, 'stage ib1': 1, 'stage ib2': 1, 'stage ic': 1, 'stage is': 1, 'is': 1,
            'stage ii': 2, 'stage iia': 2, 'stage iia1': 2, 'stage iib': 2, 'stage iic': 2, 'stage iia2': 2, 'stage iib2': 2, 'stage iic2': 2, 'i/ii nos': 2,
            'stage iii': 3, 'stage iiia': 3, 'stage iiib': 3, 'stage iiic': 3, 'stage iiic1': 3, 'stage iiic2': 3,
            'stage iv': 4, 'stage iva': 4, 'stage ivb': 4, 'stage ivc': 4,
            'stage 0': np.nan, 'stage x': np.nan, 'nan': np.nan, 'not reported': np.nan, '[not applicable]': np.nan
        }
        parsed = []
        for __ in stage_col:
            if type(__) == list:
                __ = ''.join(map(str, __))
            __ = STAGE_DICT[str(__).lower()]
            parsed.append(pd.to_numeric(__, errors='coerce'))
        
        return parsed

    def cancer_stage(self, data, impute=True, stage_cols=['clinical_stage', 'tumor_stage_diagnoses', 'ajcc_pathologic_tumor_stage']):
        """For two variables used to record stage in TCGA, parse contents and harmonize."""
        if data.filter(stage_cols).shape[1] > 1:
            print("Harmonizing multiple stage variables by taking mean.")
        sample_stages = data.filter(stage_cols).apply(self._parse_stage_col).mean(axis=1)
        if impute is True:
            sample_stages.fillna(sample_stages.mean(), inplace=True)
        data['tumor_clinical_stage'] = sample_stages

        return data.drop(stage_cols, axis=1, errors='ignore')
    
    def prepare_tcga_phenotype_data(self, data) -> pd.DataFrame:
        """Prepare data for further analysis."""
        return (data
                .pipe(self.filter_cols)
                .pipe(column_janitor)
                .rename(columns={'submitter_id_samples': 'sample'})
                .pipe(self.cancer_stage))


class MappingSampleCovariates(TCGAfromCache):
    """Prepare sample covariates data for mapping"""
    def __init__(
            self,
            cancer: str,
            covars: List[str] = ['batch_number', 'gender_demographic', 'age_at_diagnosis_diagnoses'],
            cache_dir: Path = Path("resources")
        ) -> None:
        self.cache_dir = cache_dir
        self.pheno_covars = covars
        self.ancestry_pcs = TCGAancestryPCs(cancer).ancestry_pcs
        super().__init__(cancer=cancer, cache_dir=cache_dir)
        
        self.sample_covar_cache = self.cache_dir / f'cancer~{cancer}/sample/sample_covariates_mapping.feather'
        if not self.sample_covar_cache.exists():
            print('Preparing sample covariates data ...')
            self.sample_covar_cache.parent.mkdir(parents=True, exist_ok=True)
            self.prepare_sample_covariate_data().to_feather(self.sample_covar_cache)

    @property
    def phenotype(self):
        return (self.fetch_from_omic('phenotype')
                .filter(['sample'] + self.pheno_covars)
                .dropna(how='all', axis=1))

    @property
    def shared_samples(self) -> Dict:
        return utils.fuzzymatch_samples(key_samples=self.ancestry_pcs['sample'], val_samples=self.phenotype['sample'])

    @property
    def prepared_data(self) -> pd.DataFrame:
        return self().set_index('sample')

    def prepare_sample_covariate_data(self) -> pd.DataFrame:
        return (self.ancestry_pcs
                .eval("sample = sample.map(@self.shared_samples.get)")
                .dropna()
                .set_index('sample')
                .add_prefix('ancestry_')
                .reset_index()
                .merge(self.phenotype, how='outer')
                .groupby('sample', group_keys=False)
                .apply(lambda df: df.mode(numeric_only=False, dropna=True))
                .dropna()
                .reset_index(drop=True))

    def prepare_phenotype_data_for_mapping(self):
        return self.phenotype.filter(['sample'] + self.pheno_covars).dropna(how='all', axis=1)

    def __call__(self):
        return pd.read_feather(self.sample_covar_cache)
