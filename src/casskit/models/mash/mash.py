import csv
from dataclasses import dataclass, field
import io
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Callable, Dict, List, Optional

import pandas as pd

# from .config import *
from .utils import make_mash_rdata, prepare_mash_object, stratified_sample
from casskit import config, utils
from casskit.io import utils as io_utils


# Doc links
_fastqtltomash = "https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb"
_mashr17 = "https://github.com/stephenslab/mashr/issues/17"
_ashr76 = "https://github.com/stephens999/ashr/issues/76"
_mashnullcorr1 = "https://stephenslab.github.io/mashr/articles/intro_correlations.html"
_mashnullcorr2 = "https://stephenslab.github.io/mashr/reference/estimate_null_correlation_simple.html"

# R scripts
parent_dir = Path(__file__).parent
_MASH_NULLCORR_R = parent_dir / "R" / "nullcorr.R"
_MASH_MASHDATA_R = parent_dir / "R" / "mashdata.R"
_MASH_COVARIANCE_R = parent_dir / "R" / "covariance.R"
_MASH_MIXTURE_R = parent_dir / "R" / "mixture.R"
_MASH_STRONG_R = parent_dir / "R" / "strong.R"
_MASH_SUMMARY_R = parent_dir / "R" / "summary.R"


class Mash:
    f"""Fits a MASH model to harmonized eQTL data.

    For missing shat values, replace with a large value (e.g. 1e6). In
    {_fastqtltomash}, 1e3 is used. Hard-coding _HANDLE_NAN_S = 1e3. See
    also {_mashr17} and {_ashr76}.

    Args:
      tissues:
        tissues or conditions.

    Returns:
      mashed data.

    Raises:
        IOError: An error occurred accessing
    """
    __instance = None
    
    def __new__(cls, *args, **kwargs):
        # Safety check to prevent accidental instantiation
        # Singleton assumes meta-analysis will use all available data
        # (and should therefore be instantiated only once).
        if not Mash.__instance:
            Mash.__instance = object.__new__(cls)
        return Mash.__instance

    def __init__(
        self,
        tissues,
        cache_dir: Optional[Path] = None,
        rand_sample_size: int = 2E5,
        strong_sample_size: int = 2E4,
        strong_sample_method: str = "cis",
        null_corr_method: str = "simple",
    ):
        # Caching
        self.cache_dir = cache_dir if cache_dir else Path(config.CACHE_DIR)
        
        # Mash parameters
        self.tissues = tissues
        self.rand_sample_size = rand_sample_size
        self.strong_sample_size = strong_sample_size
        self.strong_sample_method = strong_sample_method
        self.null_corr_method = null_corr_method

        self.fitted = False

    def fit(self, X, strat_idx: List = ["gene_id", "Chromosome", "arm"]) -> None:
        if self.fitted is False:
            self._check_X(X)
            self._fit(X, strat_idx)
            self.fitted = True

        return self

    def transform(self, X, strat_idx: List = ["gene_id", "Chromosome", "arm"]) -> pd.DataFrame:
        if self.fitted is False:
            raise ValueError("Mash model must be fit before transforming data")

        # Temporary files for mash
        with tempfile.TemporaryDirectory() as d:
            mash_rds = f"{d}/data.rds"
            
            # Prepare data
            data_idx = prepare_mash_object(df=X,
                                           strat_idx=strat_idx,
                                           rscript=_MASH_MASHDATA_R,
                                           temp_dir=d,
                                           v_rds=getattr(self.nullcorr, "cache"),
                                           out_rds=mash_rds,
                                           condition_cols=self.tissues,
                                           )
            
            # Run mash
            mash_output = self._transform(mash_rds, getattr(self.mixture, "cache"))

        return mash_output.set_index(data_idx)

    def _check_X(self, X):
        for col in ["tissue", "Chromosome", "cis", "gene_id", "arm", "Start", "End", "coef", "std_err"]:
            if col not in X.columns:
                raise ValueError(f"X must have column {col}")

    def _fit(self, data, strat_idx) -> None:
        with tempfile.TemporaryDirectory() as d:
            
            # (1) Temporary files for mash
            # ----------------------------
            print("Creating temporary files for mash")
            self.subset_files = MashSubsetFiles(d)
            self.subset_raw_data = \
                MashSubsetData(harmonized=data,
                               strat_idx=strat_idx,
                               strong_method=self.strong_sample_method,
                               rand_size=self.rand_sample_size,
                               strong_size=self.strong_sample_size,
                               )

            for k, v in vars(self.subset_files).items():
                if ".csv" in v:
                    getattr(self.subset_raw_data, k).to_csv(v, index=False)

            # (2) Calculate cross-tissue correlation from null tests
            # ------------------------------------------------------
            print("Calculating cross-tissue correlation from null tests")
            self.nullcorr = MashNullCorr(self.null_corr_method,
                                         self.cache_dir,
                                         )
            self.nullcorr.fit(self.subset_files)

            # (3) Set subset data objects with correlation structure
            # ------------------------------------------------------
            print("Setting subset data objects with correlation structure")
            self._make_fit_data()

            # (4) Calculate data-driven covariance matrices
            # ---------------------------------------------
            print("Calculating data-driven covariance matrices")
            self.covariance = MashDataDrivenCovariance(self.cache_dir)
            self.covariance.fit(self.subset_files)

            # (5) Fit mixture model
            # ---------------------
            print("Fitting mixture model")
            self.mixture = MashMixtureMod(getattr(self.covariance, "cache"),
                                          self.cache_dir,
                                          )
            self.mixture.fit(self.subset_files)

    def _make_fit_data(self):
        for subset in ["random", "strong"]:
            make_mash_rdata(_MASH_MASHDATA_R,
                            getattr(self.subset_files, f"bhat_{subset}"),
                            getattr(self.subset_files, f"shat_{subset}"),
                            getattr(self.nullcorr, "cache"),
                            getattr(self.subset_files, f"hat_{subset}"),
                            )

    def _transform(self, mash_data, mix_prop):
        STDOUT_MSG_LINES = 5 # Number of lines to skip in stdout
        
        # Run mash
        mash_l = []
        with subprocess.Popen(["Rscript", _MASH_SUMMARY_R,
                               "--data", mash_data,
                               "--mix", mix_prop
                               ],
                              stdout=subprocess.PIPE
                              ) as p:
            with io.TextIOWrapper(p.stdout, newline=os.linesep) as f:
                reader = csv.reader(f, delimiter=",")
                for r in reader:
                    mash_l.append(pd.Series(r))
        
        # Convert to dataframe
        return pd.concat(mash_l[STDOUT_MSG_LINES:], axis=1).set_index(0).T.dropna(how="all")

    def __repr__(self):
        return f"Mash(fitted={self.fitted})"

    def __str__(self):
        return ("mash model")

####################################################################################################
# This needs to be refactored

@dataclass
class MashSubsetFiles:
    """Paths to temporary subset files."""
    CSV_FILES = ["bhat_random", "shat_random", "bhat_strong", "shat_strong"]
    RDS_FILES = ["hat_random", "hat_strong"]
    
    temp_dir: Path

    def __post_init__(self):
        if not isinstance(self.temp_dir, Path):
            self.temp_dir = Path(self.temp_dir)

        for __ in self.CSV_FILES:
            setattr(self, __, (self.temp_dir / f"{__}.csv").as_posix())

        for __ in self.RDS_FILES:
            setattr(self, __, (self.temp_dir / f"{__}.rds").as_posix())

        self.temp_dir = self.temp_dir.as_posix()

@dataclass
class MashSubsetData:
    
    # Parameters
    _HANDLE_NAN_B = 0
    _HANDLE_NAN_S = 1E3
    _B = "coef"
    _S = "std_err"
    
    # Args
    harmonized: pd.DataFrame
    strat_idx: list[str] = field(default_factory=list)
    rand_size: int = 2E5
    strong_method: str = "cis"
    strong_size: int = 2E4

    # Cache
    read_cache: Callable = lambda cache: pd.read_csv(cache)
    write_cache: Callable = lambda data, cache: data.to_csv(cache, index=False)

    @io_utils.cache_on_disk
    def make_random_tests(self) -> pd.DataFrame:
        return stratified_sample(self.harmonized,
                                 classes=self.strat_idx,
                                 size=self.rand_size,
                                 by="tissue",
                                 )

    @io_utils.cache_on_disk
    def make_strong_tests(self) -> pd.DataFrame:
        if self.strong_method == "mash_1by1":
            raise NotImplementedError(f"{_MASH_STRONG_R} not implemented yet.")

        elif self.strong_method == "cis":
            return stratified_sample(self.harmonized.query("cis == True"),
                                     classes=self.strat_idx,
                                     size=self.strong_size,
                                     by="tissue",
                                     )

    def _make_mash_mats(self, df, suffix) -> None:
        """Make matrices for mash input."""
        setattr(self, f"bhat_{suffix}",
                self._make_tissue_mat(df, self._B, self._HANDLE_NAN_B))
        
        setattr(self, f"shat_{suffix}",
                self._make_tissue_mat(df, self._S, self._HANDLE_NAN_S))

    def _make_tissue_mat(self, df, val, fill) -> pd.DataFrame:
        """Make matrix for mash input."""
        return df.pivot_table(index=self.strat_idx,
                              columns="tissue",
                              values=val,
                              fill_value=fill
                              )

    def __post_init__(self):
        if self.strong_method not in ["mash_1by1", "cis"]:
            raise ValueError(f"Method {self.strong_method} not supported.")

        # Cache
        self.cache_dir = config.get_cache()
        self.read_cache = lambda cache: pd.read_csv(cache)
        self.write_cache = lambda data, cache: data.to_csv(cache, index=False)

        # Random subset of eQTLs
        self.path_cache = Path(
            self.cache_dir,
            f"map.harmonized.hat_random.size~{int(self.rand_size)}.csv"
        )
        self.hat_random = self.make_random_tests()
        self._make_mash_mats(self.hat_random, suffix="random")
        
        # Strong subset of eQTLs
        self.path_cache = Path(
            self.cache_dir,
            f"map.harmonized.hat_strong.method~{self.strong_method}.csv"
        )
        self.hat_strong = self.make_strong_tests()
        self._make_mash_mats(self.hat_strong, suffix="strong")

####################################################################################################

class MashNullCorr:
    def __init__(self, method: str = "simple", cache_dir: Optional[Path] = None):
        self.method = method
        self.cache_dir = cache_dir
        self.cache = self.cache_dir / f"mash.method~{self.method}.nullvar.rds"
        self.fitted = False

    def fit(self, subset_files: Dict) -> None:
        if self.fitted is False:
            self.calculate_null_corr(subset_files)
            self.fitted = True
        
        return self

    def calculate_null_corr(self, subset_files) -> None:
        f"""Estimate null cross-tissue correlations.

        Uses null tests to calculate covariances / correlations among
        tissues. A wrapper for the R script {_MASH_NULLCORR_R}.
        For more info, see {_mashnullcorr1} and {_mashnullcorr2}.

        Args:
          subset_files:
            A MashSubsetFiles instance.

        """
        args = dict()
        args["method"] = self.method
        args["out"] = self.cache
        args["b_rand"] = subset_files.__dict__.get("bhat_random")
        args["s_rand"] = subset_files.__dict__.get("shat_random")
        
        utils.subprocess_cli_rscript(_MASH_NULLCORR_R, args)
        
    def __repr__(self):
        return f"MashNullCorr(method={self.method!r}, fitted={self.fitted})"

    def __str__(self):
        return ("Mash, null correlation/covariance: "
                f"{self.method!r}")

class MashDataDrivenCovariance:
    def __init__(self,cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir
        self.cache = self.cache_dir / f"mash.ddcovariance.rds"
        self.fitted = False

    def fit(self, subset_files: Dict) -> None:
        if self.fitted is False:
            self.calculate_covariances(subset_files)
            self.fitted = True
        
        return self

    def calculate_covariances(self, subset_files) -> None:
        """Perform extreme deconvolution to estimate data-driven covariance matrices."""
        args = dict()
        args["out"] = self.cache
        args["strong"] = subset_files.__dict__.get("hat_strong")

        utils.subprocess_cli_rscript(_MASH_COVARIANCE_R, args)

class MashMixtureMod:
    def __init__(self, cov_rds: Path, cache_dir: Optional[Path] = None):
        self.cov_rds = cov_rds
        self.cache_dir = cache_dir
        self.cache = self.cache_dir / f"mash.mixturemod.rds"
        self.fitted = False

    def fit(self, subset_files: Dict) -> None:
        if self.fitted is False:
            self.fit_mixture_model(subset_files)
            self.fitted = True
        
        return self

    def fit_mixture_model(self, subset_files):
        args = dict()
        args["covariance"] = self.cov_rds
        args["out"] = self.cache
        args["random"] = subset_files.__dict__.get("hat_random")

        utils.subprocess_cli_rscript(_MASH_MIXTURE_R, args)