import logging
import sys

from . import io
from . import models as mod
from . import pipelines as pipe
from . import preprocess as pp
from .config import set_cache

from .models.latent.expression_pcs import BatchModelEPCS as ePCS

from .preprocess.expression import (
    CountThreshold,
    EdgeRCPM,
    ProteinCoding,
    RINT,
    ToCounts,
    VariationThreshold
)
from pipelines.gtex import GTEx

from io.annotate import build_ensembl_cache, get_ensembl, annotate_genes
from io.ontology.get_ontology import build_ontology_cache, get_ontology
from io.tcga.get_tcga import build_tcga_cache, get_tcga


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# https://xenabrowser.net/datapages/

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
    "io",
    "mod",
    "pipe",
    "pp",
    "set_cache",
    "ePCS",
    "CountThreshold",
    "EdgeRCPM",
    "ProteinCoding",
    "RINT",
    "ToCounts",
    "VariationThreshold",
    "GTEx",
    "build_ensembl_cache",
    "get_ensembl",
    "annotate_genes",
    "build_ontology_cache",
    "get_ontology",
    "build_tcga_cache",
    "get_tcga",
]