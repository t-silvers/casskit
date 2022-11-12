import pandas as pd

from .biogrid import get_biogrid
from .corum import get_corum
from .cosmic import get_cosmic
from .trrust import get_trrust


def build_ontology_cache():
    get_biogrid()
    get_corum()
    get_cosmic()
    get_trrust()
    
def get_ontology(resource) -> pd.DataFrame:
    """Get ontology / pathway / ... data from local cache."""
    match resource.lower():
        
        case "biogrid":
            return get_biogrid()

        case "corum":
            return get_corum()

        case "cosmic":
            return get_cosmic()
        
        case "trrust":
            return get_trrust()
                
        case _:
            raise ValueError(f"Resource {resource} not found.")