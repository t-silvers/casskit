import re

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