import re
from difflib import SequenceMatcher
from typing import Dict

import pandas as pd


def column_janitor(df: pd.DataFrame) -> pd.DataFrame:
    """Make all columns lowercase and convert special characters to underscores."""
    return df.rename(columns=str.lower).rename(columns=lambda x: janitor(x))

def janitor(a: str) -> str:
    """Make all values lowercase and convert special characters to underscores in string a.
    Strip trailing underscores from the string.
    
    Parameters
    ----------
    - a string
    
    Returns
    -------
    - a string
    """
    return re.sub("[^a-zA-Z0-9_]", "_", a.lower()).rstrip('_')

def fuzzymatch(k:list, v:list) -> Dict:
    """Fuzzy match IDs in (k)ey to IDs in (v)als
    
    Args
    -------
    key_samples: samples in dict key (<- "to replace...")
    val_samples: samples in dict values (<- "... with these")

    Returns
    -------
    A dict of matched names. A warning is given for un-matched samples.

    Example
    -------
    TCGA use case:
        Full TCGA sample IDs have 4-5 sections, separated by dashes. The first
        sections identify the donor, the last sections describe the sample. Many
        databases of TCGA analyses--particularly those with donor-level analyses--
        drop the last sections.

        This function matches TCGA names of varying completeness.

    """
    # Remove duplicates
    k_ = list(set(k))
    v_ = list(set(v))

    match_dict = {}
    for k in k_:
        for v in v_:
            match_size = SequenceMatcher(None, k, v).find_longest_match().size
            if (match_size == len(k) or match_size == len(v)):
                match_dict[k] = v

    return match_dict
