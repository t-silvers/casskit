# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

import pandas as pd

from casskit.io.utils import check_package_version
from casskit.io.ontology.biogrid import get_biogrid
from casskit.io.ontology.corum import get_corum
from casskit.io.ontology.cosmic import get_cosmic
from casskit.io.ontology.trrust import get_trrust


def build_ontology_cache(cosmic_data):
    get_biogrid()
    get_corum()
    get_cosmic(data=cosmic_data)
    get_trrust()
    
def get_ontology(resource) -> pd.DataFrame:
    """Get ontology / pathway / ... data from local cache."""
    resource = resource.lower()
    
    # Still hits syntax checkers
    # if check_package_version("python", "3.10") is True:

    #     match resource:
    #         case "biogrid":
    #             return get_biogrid()

    #         case "corum":
    #             return get_corum()

    #         case "cosmic":
    #             return get_cosmic()
            
    #         case "trrust":
    #             return get_trrust()
                    
    #         case _:
    #             raise ValueError(f"Resource {resource} not found.")
    
    # # For backwards compatibility with Python <=3.9
    # else:
    
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