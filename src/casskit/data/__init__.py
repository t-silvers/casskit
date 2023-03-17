from pkgutil import extend_path

from . import datasets
from .base import BaseOmic, BaseMultiOmic
from .multiomic import MultiOmic
from .omic import (
    CopyNumberVariation,
    MessengerRNA,
    Protein,
    SomaticMutation,
)

__all__ = [
    "datasets",
    "BaseOmic",
    "BaseMultiOmic",
    "MultiOmic",
    "CopyNumberVariation",
    "MessengerRNA",
    "Protein",
    "SomaticMutation",
]
__path__ = extend_path(__path__, __name__)
