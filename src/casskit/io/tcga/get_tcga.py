from ...io.tcga._gdcxena import build_tcga, get_gdc_tcga
from ...io.tcga._ancestry import get_ancestry_pcs
from ...io.tcga._subtype import get_subtypes
from ...io.tcga._tumorpurity import get_tumor_purity


def build_tcga_cache(cancer):
    """Build local cache of TCGA data."""
    build_tcga(cancer)
    get_ancestry_pcs(cache_only=True)
    get_subtypes(cache_only=True)
    
def get_tcga(feature, cancer=None):
    """Get TCGA data from local cache."""
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