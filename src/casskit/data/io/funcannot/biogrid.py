# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from dataclasses import dataclass, field
from pathlib import Path
import requests
from typing import Optional

import pandas as pd

from ..base import DataURLMixin
from ..config import CACHE_DIR
from ..utils import cache_on_disk, column_janitor


# TODO: Refactor biogrid to use REST API?
try:
    BIOGRID_ASSET = pd.read_csv(
        Path(__file__).parent / "assets/BIOGRID-MV-Physical-4.4.218.tab3.txt",
        sep="\t", low_memory=False
    )
except:
    BIOGRID_ASSET = pd.DataFrame()

# BIOGRID_URL = "https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.212/BIOGRID-ALL-4.4.212.tab3.zip"
# BIOGRID_URL = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip"
BIOGRID_URL = "https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip"


@dataclass
class BioGRID(DataURLMixin):
    """Fetch BioGRID interaction data."""
    cache_dir: Optional[Path] = field(init=True, default=CACHE_DIR)
    url: str = BIOGRID_URL
    organism: str = "Homo sapiens"
    
    @cache_on_disk
    def fetch(self) -> pd.DataFrame:
        if False:
            data = pd.read_csv(self.url, compression="zip", sep="\t")
        else:
            data = BIOGRID_ASSET

        return (data
                .pipe(column_janitor)
                .query("""
                       organism_name_interactor_a == @self.organism and \
                       organism_name_interactor_b == @self.organism
                       """)
                .filter([
                    "official_symbol_interactor_a", "official_symbol_interactor_b",
                    "synonyms_interactor_a", "synonyms_interactor_b",
                ]))

    def set_cache(self, cache_dir: Path) -> Path:
        self.path_cache = Path(cache_dir, f"biogrid.pkl")
        self.read_cache = lambda cache: pd.read_pickle(cache)
        self.write_cache = lambda data, cache: data.to_pickle(cache)

    @classmethod
    def get(cls, cache_only: bool = False) -> pd.DataFrame:
        data = cls().fetch()
        if cache_only is False:
            return data

    def __post_init__(self):
        self.set_cache(self.cache_dir)

get_biogrid = BioGRID.get
"""Convencience functions for loading BioGRID data."""

#######
# DEV #
#######

def _biogrid_rest_evidence():
    params = {"accesskey": "6f5b1db04594466337d179d76c147877", "format": "json"}
    r = requests.get("https://webservice.thebiogrid.org/evidence/", params=params)
    evidence = r.json()
    
    return list(evidence.keys())

def biogrid_rest_api(gene_list=None, exclude_genes="false", max=10000):
    request_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accesskey": "6f5b1db04594466337d179d76c147877",
        "format": "json",  # Return results in TAB2 format
        "max": max,
        "excludeGenes": exclude_genes,
        "searchNames": "true",  # Search against official names
        "includeInteractors": "true",  # Set to true to get any interaction involving EITHER gene, set to false to get interactions between genes
        "includeInteractorInteractions": "false",  # Set to true to get interactions between the geneList’s first order interactors
        "taxId": 9606,  # Limit to H Sapiens
        "throughputTag": "low",  # Limit to low throughput
    }
    if gene_list is not None:
        params.update({"geneList": "|".join(gene_list)})

    r = requests.get(request_url, params=params)
    interactions = r.json()
    
    # Create a hash of results by interaction identifier
    data = {}
    for interaction_id, interaction in interactions.items():
        data[interaction_id] = interaction
        # Add the interaction ID to the interaction record, so we can reference it easier
        data[interaction_id]["INTERACTION_ID"] = interaction_id

    # Load the data into a pandas dataframe
    dataset = pd.DataFrame.from_dict(data, orient="index")

    # Re-order the columns and select only the columns we want to see

    columns = [
        "INTERACTION_ID",
        "ENTREZ_GENE_A",
        "ENTREZ_GENE_B",
        "OFFICIAL_SYMBOL_A",
        "OFFICIAL_SYMBOL_B",
        "EXPERIMENTAL_SYSTEM",
        "PUBMED_ID",
        "PUBMED_AUTHOR",
        "THROUGHPUT",
        "QUALIFICATIONS",
    ]
    dataset = dataset[columns]

    return dataset

def build_biogrid_from_rest():
    """
    Silly double-counting here because interactions are asymmetric (ie lists A-B separate from B-A)
    """
    genes_queried = []
    
    # Get sample of biogrid interactions
    biogrid_iter0_A = biogrid_rest_api()
    gene_ids = biogrid_iter0_A["OFFICIAL_SYMBOL_A"].unique().tolist()
    
    # Use sample to query
    biogrid_iter0_B = biogrid_rest_api(gene_list=gene_ids)
    gene_ids = biogrid_iter0_B["OFFICIAL_SYMBOL_A"].unique().tolist()
    genes_queried = genes_queried + gene_ids

    # Repeat
    biogrid_iter1_A = biogrid_rest_api(gene_list=genes_queried, exclude_genes="true")
    gene_ids = biogrid_iter1_A["OFFICIAL_SYMBOL_A"].unique().tolist()

    biogrid_iter1_A = biogrid_rest_api(gene_list=gene_ids)
    gene_ids = biogrid_iter1_A["OFFICIAL_SYMBOL_A"].unique().tolist()