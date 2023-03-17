from pkgutil import extend_path

from .gdc_xena import build_tcga, get_gdc_tcga
from .tcgabiolinks_subtype import get_subtypes


__all__ = ["build_tcga_cache", "get_tcga"]
__path__ = extend_path(__path__, __name__)

def build_tcga_cache(cancer, minimal=False, n_samples=20):
    """Build local cache of TCGA data."""

    build_tcga(cancer, minimal=minimal, n_samples=n_samples)
    get_subtypes(cache_only=True)

def get_tcga(data_name, cancer=None):
    """Get TCGA data from local cache."""
    
    if data_name == "subtypes":
        return get_subtypes()
    
    else:
        try:
            return get_gdc_tcga(cancer, data_name)
        except Exception as e:
            raise ValueError(f"data_name {data_name} not found.") from e