from .ontology.get import build_ontology_cache, get_ontology
from .tcga.get import build_tcga_cache, get_tcga

__path__ = __import__('pkgutil').extend_path(__path__, __name__)
__all__ = [
    "build_ontology_cache",
    "get_ontology",
    "build_tcga_cache",
    "get_tcga"
]