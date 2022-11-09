import logging
import os

os.environ["CASSKIT_CACHE_DIR"] = "/scratch/users/tsilvers/ceqtl_selection"

from ..io._base import (
    DataURLMixin,
    ElsevierLink,
)

from ..io.ontology.get_ont import get_ont
from ..io.tcga.get_tcga import build_tcga_cache, get_tcga

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# https://xenabrowser.net/datapages/

__all__ = [
    "DataURLMixin",
    "TCGAXenaLoader",
    "ElsevierLink",
    "get_ont",
    "build_tcga_cache",
    "get_tcga",
]