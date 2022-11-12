# from .annotate import build_ensembl_cache, get_ensembl, annotate_genes
# from .ontology.get_ontology import build_ontology_cache, get_ontology
# from .tcga.get_tcga import build_tcga_cache, get_tcga

# __all__ = [
#     "build_ensembl_cache",
#     "get_ensembl",
#     "annotate_genes",
#     "build_ontology_cache",
#     "get_ontology",
#     "build_tcga_cache",
#     "get_tcga",
# ]

from . import ontology
from . import tcga

__path__ = __import__('pkgutil').extend_path(__path__, __name__)

__all__ = [
    "ontology",
    "tcga"
]
