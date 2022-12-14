import logging
import sys

from .config import set_cache, set_logging
from .dask_cluster import DaskCluster
from . import data, io, factory
from . import models as mod
from . import pipelines as pipe
from . import preprocess as pp
from . import selectors as sel

dask_cluster = DaskCluster.dask_cluster

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

__path__ = __import__('pkgutil').extend_path(__path__, __name__)

__all__ = [
    "set_cache",
    "set_logging",
    "dask_cluster",
    "data",
    "io",
    "factory",
    "mod",
    "pipe",
    "pp",
    "sel",
]
