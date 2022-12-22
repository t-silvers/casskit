import csv
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import List, Optional

import casskit as ck
import numpy as np
import pandas as pd

from cneqtl.data import HarmonizedCNeQTLs


# TODO: Refactor into casskit


# Doc links
_fastqtl_to_mash = "https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb"
_mashr_issues_17 = "https://github.com/stephenslab/mashr/issues/17"
_ashr_issues_76 = "https://github.com/stephens999/ashr/issues/76"

# R scripts
if sys.stdin.isatty():
    parent_dir = Path("").resolve()

else:
    parent_dir = Path(__file__).parent

_MASH_NULLCORR_R = parent_dir / "nullcorr.R"
_MASH_DATA_R = parent_dir / "mashdata.R"
_MASH_COVARIANCE_R = parent_dir / "covariance.R"
_MASH_MIXTURE_R = parent_dir / "mixture.R"
_MASH_STRONG_R = parent_dir / "strong.R"
_MASH_SUMMARY_R = parent_dir / "summary.R"


def stratified_sample(
    df: pd.DataFrame,
    classes: List,
    size: int,
    by: str,
    samples: Optional[int] = 1E4,
    min_strat: Optional[int] = 5,
    n_strat_samples: Optional[int] = 10,
    oversample: bool = False,
) -> pd.DataFrame:
    """Stratified sample of a dataframe.
    
    Performs a version of stratified sampling useful for
    mash analysis. Sampling enforces each tissues is
    represented, without oversampling or bias for large-
    effect size cn-eQTLs.
    
    The design here is a form of adaptive, stratified,
    two-phase sampling (Thompson, 2012):
    
    1. Sample cn-eQTL-eGene pairs (stage 1)
    2. Evaluate representation of classes (stage 1)
    3. Continue sampling until all classes are represented
    4. Return a sample from stage 1 (stage 2)

    """
    if oversample is True:
        raise NotImplementedError("Oversampling not implemented.")

    # Subsampling for stratified sampling
    n_by = df[by].nunique()
    n_per_strat = df.groupby(classes).count()[by].mean()
    n_strat = int(size/(n_by/n_per_strat))

    # Stratified sampling
    # Strategy: Randomly sample n times and evaluate if classes are represented
    df_ = df.set_index(classes)
    samples_dfs = []
    samples_drawn = 0
    while len(samples_dfs) < n_strat_samples:
        # Random sampling
        idx_sample = df_.index.to_frame(index=False).sample((n_strat)).set_index(classes).index
        df_sample = df_.loc[idx_sample]
        
        # Check if classes are represented
        n_sample_classes = df_sample.tissue.value_counts().mask(lambda x: x < min_strat).dropna().size
        if n_sample_classes == n_by:
            # print("Found")
            samples_dfs.append(
                df_sample.reset_index()
            )
        
        samples_drawn += 1
        # print(samples_drawn)
        if samples_drawn > samples:
            break
    
    if len(samples_dfs) == 0:
        raise ValueError("No stratified sample could be drawn. Consider increasing "
                         "the number of samples or decreasing the minimum number "
                         "of samples per class.")
    
    rng = np.random.default_rng()
    
    return samples_dfs[rng.choice(range(len(samples_dfs)))]

class Mash:
    f"""
    Notes
    - For missing shat values, replace with a large value (e.g. 1e6). In
    {_fastqtl_to_mash}, 1e3 is used. Hard-coding _HANDLE_NAN_S = 1e3. See
    also {_mashr_issues_17} and {_ashr_issues_76}.
    """
    
    # prepare mash input
    _STRAT_CLASSES = ["gene_id", "Chromosome", "arm"]
    _HANDLE_NAN_B = 0
    _HANDLE_NAN_S = 1E3
    _B = "coef"
    _S = "std_err"
    
    def __init__(
        self,
        harmonized: HarmonizedCNeQTLs,
        rand_size: int = 2E5,
        strong_method: str = "cis",
        strong_size: int = 2E4,
        nullvar_method: str = "simple",
    ) -> None:
        self.fitted = False
        
        self.harmonized_data = harmonized.harmonized_data
        self.rand_size = rand_size
        self.strong_method = strong_method
        self.strong_size = strong_size
        self.nullvar_method = nullvar_method
        
        # Cache
        self.cache_dir = ck.config.get_cache()
        self.read_cache = lambda cache: pd.read_csv(cache)
        self.write_cache = lambda data, cache: data.to_csv(cache, index=False)
        
        # mash intermediates, save to re-use
        self.nullvar_cache, self.ddcovariance_cache = \
            [self.cache_dir / f"mash.method~{self.nullvar_method}.{__}.rds" for 
             __ in ["nullvar", "ddcovariance"]]

    @property
    def hat_random(self) -> pd.DataFrame:
        self.path_cache = Path(
            self.cache_dir,
            f"map.harmonized.hat_random.size~{int(self.rand_size)}.csv"
        )
        hat_random_tidy = self.make_random_tests()
        
        # Make data for mash
        self._make_mash_mats(hat_random_tidy, suffix="random")
        
        return hat_random_tidy

    @property
    def hat_strong(self) -> pd.DataFrame:
        self.path_cache = Path(
            self.cache_dir,
            f"map.harmonized.hat_strong.method~{self.strong_method}.csv"
        )
        hat_strong_tidy = self.make_strong_tests(self.harmonized_data,
                                                 size=self.strong_size,
                                                 method=self.strong_method,
                                                 )

        # Make data for mash
        self._make_mash_mats(hat_strong_tidy, suffix="strong")

        return hat_strong_tidy

    def _make_mash_mats(self, df: pd.DataFrame, suffix: str) -> None:
        """Make matrices for mash input."""
        setattr(self, f"bhat_{suffix}",
                self._make_tissue_mat(df, self._B, self._HANDLE_NAN_B))
        setattr(self, f"shat_{suffix}",
                self._make_tissue_mat(df, self._S, self._HANDLE_NAN_S))

    def _make_tissue_mat(self, df: pd.DataFrame, val: str, fill: str) -> pd.DataFrame:
        """Make matrix for mash input."""
        return df.pivot_table(index=self._STRAT_CLASSES,
                              columns="tissue",
                              values=val,
                              fill_value=fill
                              )

    @ck.io.utils.cache_on_disk
    def make_random_tests(self) -> pd.DataFrame:
        return stratified_sample(self.harmonized_data,
                                 classes=self._STRAT_CLASSES,
                                 size=self.rand_size,
                                 by="tissue",
                                 )

    @ck.io.utils.cache_on_disk
    def make_strong_tests(
        self,
        data: pd.DataFrame,
        size: int = 1E8,
        method: str = "cis",
    ) -> pd.DataFrame:
        if method not in ["mash_1by1", "cis"]:
            raise ValueError(f"Method {method} not supported.")

        if method == "mash_1by1":
            strong_tests_l = []
            with subprocess.Popen(["Rscript", _MASH_STRONG_R],
                                  stdout=subprocess.PIPE,
                                  ) as p:
                with io.TextIOWrapper(p.stdout, newline=os.linesep) as f:
                    reader = csv.reader(f, delimiter=",")
                    
                    for r in reader:
                        strong_tests_l.append(pd.Series(r))

            return pd.concat(strong_tests_l[2:], axis=1).set_index(0).T

        elif method == "cis":
            # return stratified_sample(data.query("cis == True"),
            #                          classes=self._STRAT_CLASSES,
            #                          size=size,
            #                          by="tissue",
            #                          )
            
            return data.query("cis == True")

    def fit(self):
        """Fit mashr model."""
        if not self.nullvar_cache.exists() or not self.ddcovariance_cache.exists():
            # Temporary subsample files for mash
            with tempfile.TemporaryDirectory() as d:
                temp_dir = Path(self.cache_dir, d)
                temp_fns = dict(bhat_random=temp_dir / "bhat.random.csv",
                                shat_random=temp_dir / "shat.random.csv",
                                hat_random=temp_dir / "random.rds",
                                bhat_strong=temp_dir / "bhat.strong.csv",
                                shat_strong=temp_dir / "shat.strong.csv",
                                hat_strong=temp_dir / "strong.rds",
                                )

                # Temporary files for mash
                for k, v in temp_fns.items():
                    getattr(self, k).to_csv(v, index=False)

                # (1) Calculate cross-tissue correlation from null tests
                self.calculate_null_corr(temp_fns["bhat_random"],
                                         temp_fns["bhat_strong"],
                                         method="simple"
                                         )

                # (2) Set data objects with correlation structure
                self.set_data_objects(temp_fns)

                # (3) Calculate data-driven covariance matrices
                self.calculate_covariances()
                
                # (4) Fit mixture model
                self.fit_mixture_model()

            # Remove temporary directory
            shutil.rmtree(d)
        
        self.fitted = True
        
        # (5) Calculate mash results
        self.mash_summary()

        return self

    def transform(self, harmonized: HarmonizedCNeQTLs):
        """Run mashr."""
        if not self.fitted:
            raise ValueError("Model not fitted.")

    def calculate_null_corr(self, fn_dict) -> None:
        cmd = (f"Rscript {_MASH_NULLCORR_R} "
               f" --method {self.nullvar_method} --out {self.nullvar_cache} "
               "--b_rand {bhat_random} --s_rand {shat_random}".format(**fn_dict))

        subprocess.run(cmd, shell=True, check=True)

    def set_data_objects(self, fn_dict) -> None:
        cmd = (f"Rscript {_MASH_DATA_R} --nullcorr {self.nullvar_cache} "
               "--b_rand {bhat_random} --s_rand {shat_random} "
               "--out_rand {hat_random} --b_strong {bhat_strong} "
               "--s_strong {shat_strong} --out_strong {hat_strong}"
               .format(**fn_dict))

        subprocess.run(cmd, shell=True, check=True)
    
    def calculate_covariances(self, fn_dict) -> None:
        """Perform extreme deconvolution to estimate data-driven covariance matrices."""
        if not self.nullvar_cache:
            raise ValueError("Null correlation matrix not calculated.")
        
        cmd = (f"Rscript {_MASH_COVARIANCE_R} --out {self.ddcovariance_cache}"
                "--strong {hat_strong}".format(**fn_dict))
        
        subprocess.run(cmd, shell=True, check=True)

    def fit_mixture_model(self, fn_dict):
        if not self.ddcovariance_cache:
            raise ValueError("Data-driven covariance matrix not calculated.")

        cmd = (f"Rscript {_MASH_MIXTURE_R} --out {self.ddcovariance_cache}"
                "--strong {hat_strong}".format(**fn_dict))
        
        subprocess.run(cmd, shell=True, check=True)

    def mash_summary(self):
        subtypes_l = []
        with subprocess.Popen(
            ["Rscript", _MASH_SUMMARY_R],
            stdout=subprocess.PIPE
        ) as p:
            with io.TextIOWrapper(p.stdout, newline=os.linesep) as f:
                reader = csv.reader(f, delimiter=",")
                
                for r in reader:
                    subtypes_l.append(pd.Series(r))

        return pd.concat(subtypes_l[2:], axis=1).set_index(0).T
    