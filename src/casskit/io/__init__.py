import logging

from casskit.io._base import (
    DataURLMixin,
    ElsevierLink,
)

from casskit.io import get_tcga

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# https://xenabrowser.net/datapages/

__all__ = [
    "DataURLMixin",
    "TCGAXenaLoader",
    "ElsevierLink",
    "get_tcga",
]