import casskit
from casskit.io import (
    get_tumor_purity,
    get_ancestry_pcs,
    get_gdc_tcga,
)

def get_tcga(cancer, feature, cache_dir=casskit.CACHE_DIR):
    """Get TCGA data from local cache."""
    get_gdc_tcga(cancer, feature)
    get_tumor_purity()
    get_ancestry_pcs()