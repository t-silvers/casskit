from pathlib import Path
import sys


# Doc links
_fastqtltomash = "https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb"
_mashr17 = "https://github.com/stephenslab/mashr/issues/17"
_ashr76 = "https://github.com/stephens999/ashr/issues/76"
_mashnullcorr1 = "https://stephenslab.github.io/mashr/articles/intro_correlations.html"
_mashnullcorr2 = "https://stephenslab.github.io/mashr/reference/estimate_null_correlation_simple.html"

# R scripts
if sys.stdin.isatty():
    parent_dir = Path("").resolve()

else:
    parent_dir = Path(__file__).parent

_MASH_NULLCORR_R = parent_dir / "nullcorr.R"
_MASH_MASHDATA_R = parent_dir / "mashdata.R"
_MASH_COVARIANCE_R = parent_dir / "covariance.R"
_MASH_MIXTURE_R = parent_dir / "mixture.R"
_MASH_STRONG_R = parent_dir / "strong.R"
_MASH_SUMMARY_R = parent_dir / "summary.R"
