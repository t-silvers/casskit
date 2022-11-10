import sys

from .io.ontology.get_ont import get_ont
from .io.tcga.get_tcga import build_tcga_cache, get_tcga
from .struct._tiledb import DataDB

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# https://xenabrowser.net/datapages/

from .config import *

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
    "DataDB",
    "get_ont",
    "build_tcga_cache",
    "get_tcga",
]