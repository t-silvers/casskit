import logging
from pkgutil import extend_path
import sys

from . import pp, data

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

__all__ = ["pp", "data"]
__path__ = extend_path(__path__, __name__)

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

