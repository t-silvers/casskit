# from .annotate import build_ensembl_cache, get_ensembl, annotate_genes
# from .ontology.get_ontology import build_ontology_cache, get_ontology
# from .tcga.get_tcga import build_tcga_cache, get_tcga

from . import annotate
from .ontology import get_ontology
from .tcga import get_tcga

build_ensembl_cache = annotate.build_ensembl_cache
get_ensembl = annotate.get_ensembl
annotate_genes = annotate.annotate_genes

build_ontology_cache = get_ontology.build_ontology_cache
get_ontology = get_ontology.get_ontology

build_tcga_cache = get_tcga.build_tcga_cache
get_tcga = get_tcga.get_tcga


__all__ = [
    "build_ensembl_cache",
    "get_ensembl",
    "annotate_genes",
    "build_ontology_cache",
    "get_ontology",
    "build_tcga_cache",
    "get_tcga",
]
