import os
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

HOME = os.environ["HOME"]

try:
    SCRATCH = os.environ["SCRATCH"]
except KeyError:
    pass

_LOG_DIR_ENV_KEY = "CASSKIT_LOG_DIR"

try:
    LOG_DIR = Path(os.environ[_LOG_DIR_ENV_KEY]) / ".casskit"
except:
    logger.warning(f"Environment variable {_LOG_DIR_ENV_KEY} not set. Using default cache directory.")
    LOG_DIR = Path("~/.cache").expanduser()

LOG_DIR.mkdir(parents=True, exist_ok=True)

def set_logging(logging_dir, subdirs=[".casskit"]):
    """Set logging directory."""
    global LOG_DIR
    LOG_DIR = Path(logging_dir)
    if len(subdirs) > 0:
        LOG_DIR = Path(logging_dir, *subdirs)
    
    LOG_DIR.mkdir(parents=True, exist_ok=True)

def get_logging():
    """Get logging directory."""
    return LOG_DIR


_CACHE_BASE_SUBDIR = "casskit"
_CACHE_DIR_ENV_KEY = "CASSKIT_CACHE_DIR"

try:
    CACHE_DIR = Path(os.environ[_CACHE_DIR_ENV_KEY]) / ".cache" / _CACHE_BASE_SUBDIR
except:
    logger.warning(f"Environment variable {_CACHE_DIR_ENV_KEY} not set. Using default cache directory.")
    CACHE_DIR = Path("~/.cache").expanduser() / _CACHE_BASE_SUBDIR

CACHE_DIR.mkdir(parents=True, exist_ok=True)


def set_cache(cache_dir, subdirs=[".cache", _CACHE_BASE_SUBDIR]):
    """Set cache directory."""
    global CACHE_DIR
    CACHE_DIR = Path(cache_dir)
    if len(subdirs) > 0:
        CACHE_DIR = Path(cache_dir, *subdirs)
    
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    
def get_cache():
    """Get cache directory."""
    return CACHE_DIR

def set_threads(
    threads: int = None, env_var: str = "SLURM_CPUS_PER_TASK"
) -> None:
    
    if threads is None:
        try:
            threads = os.environ[env_var]
        except:
            threads = "1"

    os.environ["BLIS_NUM_THREADS"] = threads
    os.environ["MKL_NUM_THREADS"] = threads
    os.environ["NUMEXPR_NUM_THREADS"] = threads
    os.environ["OMP_NUM_THREADS"] = threads
    os.environ["OPENBLAS_NUM_THREADS"] = threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = threads
