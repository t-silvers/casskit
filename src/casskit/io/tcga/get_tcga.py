from ...io.tcga._gdcxena import build_tcga, get_gdc_tcga
from ...io.tcga._ancestry import get_ancestry_pcs
from ...io.tcga._tumorpurity import get_tumor_purity


def build_tcga_cache(cancer):
    """Build local cache of TCGA data."""
    build_tcga(cancer)
    get_ancestry_pcs(cache_only=True)
    
def get_tcga(cancer, feature):
    """Get TCGA data from local cache."""
    get_gdc_tcga(cancer, feature)
    get_tumor_purity()
    get_ancestry_pcs()