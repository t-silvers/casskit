from .gdc_xena import build_tcga, get_gdc_tcga
from .ancestry import get_ancestry_pcs
from .subtype import get_subtypes
from .tumorpurity import get_tumor_purity

# from casskit.io.tcga.gdc_xena import build_tcga, get_gdc_tcga
# from casskit.io.tcga.ancestry import get_ancestry_pcs
# from casskit.io.tcga.subtype import get_subtypes
# from casskit.io.tcga.tumorpurity import get_tumor_purity


def build_tcga_cache(cancer):
    """Build local cache of TCGA data."""

    build_tcga(cancer)
    get_ancestry_pcs(cache_only=True)
    get_subtypes(cache_only=True)

    
def get_tcga(feature, cancer=None):
    """Get TCGA data from local cache."""
    # TODO: filter samples by cancer type
    match feature:
        
        case "ancestry":
            return get_ancestry_pcs()
        
        case "purity":
            return get_tumor_purity()
        
        case "subtypes":
            return get_subtypes()
        
        case _:
            try:
                return get_gdc_tcga(cancer, feature)
            except Exception as e:
                raise ValueError(f"Feature {feature} not found.") from e