import re
from difflib import SequenceMatcher
from typing import Dict, List

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

def fuzzy_match(k: List, v: List, deduped: bool = False) -> Dict:
    """Fuzzy match IDs in (k)ey to IDs in (v)als
        
    Args
    -------
    k : List
        Samples in Dict keys (<- "to replace...")
    v : List
        Samples in Dict values (<- "... with these")
    deduped : boolean, default=False
        Is input unique?

    Returns
    -------
    A dict of matched names. A warning is given for un-matched samples.

    Notes
    -----
    Lightweight, fast version of fuzzy matching.

    For a more robust version, see:
    [dirty_cat.fuzzy_join]
    (https://github.com/dirty-cat/dirty_cat/blob/a920c4761c5f0978056c528176895c747ad6f713/dirty_cat/_fuzzy_join.py#L24)

    TCGA use case:

    Full TCGA sample IDs have 4-5 sections, separated by dashes. The first
    sections identify the donor, the last sections describe the sample. Many
    databases of TCGA analyses--particularly those with donor-level analyses--
    drop the last sections.

    This function matches TCGA names of varying completeness.

    Example
    -------

    """
    k_, v_ = k, v

    if deduped is False:
        k_ = list(set(k))
        v_ = list(set(v))

    match_dict = {}
    for k in k_:
        for v in v_:
            match_size = (SequenceMatcher(None, k, v)
                          .find_longest_match()
                          .size)
            
            if (match_size == len(k) or match_size == len(v)):
                match_dict[k] = v

    return match_dict