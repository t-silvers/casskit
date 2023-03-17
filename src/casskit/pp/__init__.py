from pkgutil import extend_path

from . import copynumber, expression
from .expression import GTEx

__all__ = ["copynumber", "expression", "GTEx"]
__path__ = extend_path(__path__, __name__)