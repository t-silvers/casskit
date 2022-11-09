from ...io.tcga._gdcxena import get_gdc_tcga
from ...io.tcga._ancestry import get_ancestry_pcs
from ...io.tcga._tumorpurity import get_tumor_purity
from ...config import CACHE_DIR # TEMP


def get_tcga(cancer, feature, cache_dir=CACHE_DIR):
    """Get TCGA data from local cache."""
    get_gdc_tcga(cancer, feature)
    get_tumor_purity()
    get_ancestry_pcs()