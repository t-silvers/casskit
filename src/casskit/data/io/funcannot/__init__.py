# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT
from pkgutil import extend_path

import pandas as pd

from .biogrid import get_biogrid
from .corum import get_corum
from .cosmic import get_cosmic
from .trrust import get_trrust


__all__ = ["build_funcannot_cache", "get_funcannot"]
__path__ = extend_path(__path__, __name__)


def build_funcannot_cache():
    get_biogrid(cache_only=True)
    get_corum(cache_only=True)
    get_cosmic(cache_only=True)
    get_trrust(cache_only=True)

def get_funcannot(resource) -> pd.DataFrame:
    """Get funcannot / pathway / ... data from local cache."""
    resource = resource.lower()

    if resource == "biogrid":
        return get_biogrid()
    
    elif resource == "corum":
        return get_corum()
    
    elif resource == "cosmic":
        return get_cosmic()
    
    elif resource == "trrust":
        return get_trrust()
    
    else:
        raise ValueError(f"Resource {resource} not found.")