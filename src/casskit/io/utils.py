from difflib import SequenceMatcher
from functools import wraps
from pathlib import Path
import pkg_resources
import platform
import re
from typing import Callable, Dict, List
import warnings

import pandas as pd
import tqdm


def cache_on_disk(f: Callable) -> Callable:
    """Cache function output on disk.
    
    Light-weight caching decorator that caches class function output on disk.
    This decorator will only work on functions that return a single pandas
    DataFrame object.
    """
    @wraps(f)
    def wrapper(self, *args, **kwargs):
        if hasattr(self, "path_cache"):
            cache = self.path_cache

            if (Path(cache).exists() and Path(cache).stat().st_size > 0):
                print(f"Loading from cache: {cache}")
                data = self.read_cache(cache)
                
            else:
                print(f"Caching to disk: {cache}")
                Path(cache).parent.mkdir(exist_ok=True, parents=True)
                data = f(self, *args, **kwargs)
                if (data is not None and data.empty is False):
                    self.write_cache(data, cache)
                else:
                    warnings.warn(f"Data is empty. Not writing to cache: {cache}")

        else:
            data = f(self, *args, **kwargs)

        return data
    return wrapper

def check_package_version(package: str, version: str = None) -> bool:
    """Check if package is installed and at least a certain version."""
    try:
        if package == "python":
            if version is None:
                return True
            else:
                py_ver = pkg_resources.parse_version(
                    platform.python_version()
                )
                return py_ver >= pkg_resources.parse_version(version)
        
        pkg_distrib = pkg_resources.get_distribution(package)
        if version is not None:
            return pkg_resources.parse_version(pkg_distrib.version) >= pkg_resources.parse_version(version)
    
    except pkg_resources.DistributionNotFound:
        return False

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
        Samples in Dict values (<- "... with these").
        Should be the shorter list.
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
        k_, v_ = list(set(k)), list(set(v))

    # Remove nan
    rm_nan = lambda x: [i for i in x if pd.notnull(i)]
    k_, v_ = rm_nan(k_), rm_nan(v_)
        
    match_dict = {}
    for v in tqdm.tqdm(v_):
        for k in k_:
            seq_match = SequenceMatcher(None, k, v)
            match_size = seq_match.find_longest_match().size
            
            if (match_size == len(k) or match_size == len(v)):
                match_dict[k] = v
                break

    return match_dict