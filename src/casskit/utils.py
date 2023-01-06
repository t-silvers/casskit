import csv
from difflib import SequenceMatcher
import io
import itertools
import os
from pathlib import Path
import re
import subprocess
from typing import Dict, List, Union

import pandas as pd
import tqdm


def pipe_concat(*args):
    return pd.concat(args, axis=0)

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

def subprocess_cli_rscript(
    script: Path,
    args: Dict,
    ret: bool = False,
    skip_lines: int = 0,
) -> Union[None, pd.DataFrame]:
    """Run a command line interface (CLI) command in a subprocess.

    Parameters
    ----------
    cmd : str
        A command line interface (CLI) command.

    Returns
    -------
    None
    """
    cmd = f"Rscript {script} "

    for k, v in args.items():
        cmd += f"--{k} {v} "

    if ret is False:
        subprocess.run(cmd, shell=True, check=True)
    
    else:
        stdout = []
        with subprocess.Popen(
            cmd, stdout=subprocess.PIPE, shell=True
        ) as proc:
            with io.TextIOWrapper(proc.stdout, newline=os.linesep) as f:
                
                # Must write stdout lines as csv, eg w/ 
                # writeLines(readr::format_csv(mash.res), stdout())
                reader = csv.reader(f, delimiter=",")
                for r in reader:
                    stdout.append(pd.Series(r))
        
        # Convert to dataframe
        return pd.concat(stdout[skip_lines:], axis=1).set_index(0).T
