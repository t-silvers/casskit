import logging

from casskit.io._base import (
    DataURLMixin,
    ElsevierLink,
)

from casskit.io.tcga._gdc_via_xena import (
    TCGAXenaLoader,
    build_tcga,
    get_tcga
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# https://xenabrowser.net/datapages/

__all__ = [
    "DataURLMixin",
    "TCGAXenaLoader",
    "ElsevierLink",
    "build_tcga",
    "get_tcga"
]