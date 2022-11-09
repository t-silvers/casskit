import os
from pathlib import Path
import logging


logger = logging.getLogger(__name__)


HOME = os.environ["HOME"]


CACHE_BASE_SUBDIR = "casskit"
CACHE_DIR_ENV_KEY = "CASSKIT_CACHE_DIR"

try:
    CACHE_DIR = Path(os.environ[CACHE_DIR_ENV_KEY]) / ".cache" / CACHE_BASE_SUBDIR
except:
    logger.warning(f"Environment variable {CACHE_DIR_ENV_KEY} not set. Using default cache directory.")
    CACHE_DIR = Path("~/.cache").expanduser() / ".cache" / CACHE_BASE_SUBDIR

CACHE_DIR.mkdir(parents=True, exist_ok=True)


UNITS = {
    "log1p": "log1p",
    "log2": "log2",
    "log10": "log10",
    "log2(count+1)": "log1p",
    "absolute": "counts",
    "abs": "counts",
    "counts": "counts"
}