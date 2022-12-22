from pathlib import Path
import subprocess
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def prepare_mash_object(
    df: pd.DataFrame,
    strat_idx: List,
    rscript: str,
    temp_dir: Path,
    v_rds: str,
    out_rds: str,
    condition_col: str = "tissue",
    param_cols: Dict[str, str] = dict(b="coef", sd="std_err"),
    fillna: Dict[str, float] = dict(b=0, sd=1E3),
) -> pd.Index:
    
    # Convert tidy data to separate matrices
    for k, param in param_cols.items():
        df_mat = df.pivot_table(index=strat_idx,
                                columns=condition_col,
                                values=param,
                                fill_value=fillna.get(k, None),
                                )
        
        df_mat.to_csv(f"{temp_dir}/{k}.csv", index=False)
    
    # Make mash data object
    make_mash_rdata(rscript, f"{temp_dir}/b.csv", f"{temp_dir}/sd.csv", v_rds, out_rds)

    # Return index of cn-eQTLs
    return df_mat.index

def make_mash_rdata(
    rscript: str,
    b_data: str,
    s_data: str,
    v_rds: str,
    out_rds: str,
) -> None:
    cmd = (f"Rscript {rscript} "
           f"--nullcorr {v_rds} "
           f"--b {b_data} "
           f"--s {s_data} "
           f"--out {out_rds} "
           )

    subprocess.run(cmd, shell=True, check=True)

def stratified_sample(
    df: pd.DataFrame,
    classes: List,
    size: int,
    by: str,
    samples: Optional[int] = 1E4,
    min_strat: Optional[int] = 5,
    n_strat_samples: Optional[int] = 10,
    oversample: bool = False,
) -> pd.DataFrame:
    """Stratified sample of a dataframe.
    
    Performs a version of stratified sampling useful for
    mash analysis. Sampling enforces each tissues is
    represented, without oversampling or bias for large-
    effect size cn-eQTLs.
    
    The design here is a form of adaptive, stratified,
    two-phase sampling (Thompson, 2012):
    
    1. Sample cn-eQTL-eGene pairs (stage 1)
    2. Evaluate representation of classes (stage 1)
    3. Continue sampling until all classes are represented
    4. Return a sample from stage 1 (stage 2)

    """
    if oversample is True:
        raise NotImplementedError("Oversampling not implemented.")

    # Subsampling for stratified sampling
    n_by = df[by].nunique()
    n_per_strat = df.groupby(classes).count()[by].mean()
    n_strat = int(size/(n_by/n_per_strat))

    # Stratified sampling
    # Strategy: Randomly sample n times and evaluate if classes are represented
    df_ = df.set_index(classes)
    samples_dfs = []
    samples_drawn = 0
    while len(samples_dfs) < n_strat_samples:
        # Random sampling
        idx_sample = df_.index.to_frame(index=False).sample((n_strat)).set_index(classes).index
        df_sample = df_.loc[idx_sample]
        
        # Check if classes are represented
        n_sample_classes = df_sample.tissue.value_counts().mask(lambda x: x < min_strat).dropna().size
        if n_sample_classes == n_by:
            # print("Found")
            samples_dfs.append(
                df_sample.reset_index()
            )
        
        samples_drawn += 1
        # print(samples_drawn)
        if samples_drawn > samples:
            break
    
    if len(samples_dfs) == 0:
        raise ValueError("No stratified sample could be drawn. Consider increasing "
                         "the number of samples or decreasing the minimum number "
                         "of samples per class.")
    
    rng = np.random.default_rng()
    
    return samples_dfs[rng.choice(range(len(samples_dfs)))]
