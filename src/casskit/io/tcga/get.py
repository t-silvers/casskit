# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from casskit.io.tcga.ancestry import get_ancestry_pcs
from casskit.io.tcga.aneuploidy_score import get_tcga_aneuploidy_scores
from casskit.io.tcga.gdc_xena import build_tcga, get_gdc_tcga
from casskit.io.tcga.subtype import get_subtypes
from casskit.io.tcga.survival_cdr import get_tcga_cdr_survival
from casskit.io.tcga.tumorpurity import get_tumor_purity


def build_tcga_cache(cancer):
    """Build local cache of TCGA data."""

    build_tcga(cancer)
    get_ancestry_pcs(cache_only=True)
    __ = get_tcga_aneuploidy_scores()
    __ = get_tumor_purity()
    get_subtypes(cache_only=True)
    __ = get_tcga_cdr_survival()
    
    
def get_tcga(feature, cancer=None):
    """Get TCGA data from local cache."""
    # TODO: filter samples by cancer type
    match feature:
        
        case "ancestry":
            return get_ancestry_pcs()

        case "aneuploidy_score":
            return get_tcga_aneuploidy_scores()
        
        case "purity":
            return get_tumor_purity()
        
        case "subtypes":
            return get_subtypes()

        case "survival":
            return get_tcga_cdr_survival()
        
        case _:
            try:
                return get_gdc_tcga(cancer, feature)
            except Exception as e:
                raise ValueError(f"Feature {feature} not found.") from e