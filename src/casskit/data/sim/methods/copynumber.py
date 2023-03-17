import numpy as np
import pandas as pd
import tqdm

from .config import CHROMOSOME_SIZES
from .generics import simulate_matrix, take_subset
from ..utils import simulate_sample_ids


def chromswap_augmented(
    n_samples: int,
    *args,
    template_copynumber: pd.DataFrame = None,
    random_state: int = None,
    **kwargs,
):
    rng = np.random.default_rng(random_state)
    req_cols = ["sample_id", "Chromosome"]
    if not all([col in template_copynumber.columns for col in req_cols]):
        raise ValueError(f"template_copynumber must have columns {req_cols}")

    # Possible values
    sample_ids = template_copynumber["sample_id"].unique()
    chromosomes = template_copynumber["Chromosome"].unique()

    # Randomly sample sample_ids with replacement
    samples_swapped_np = rng.choice(sample_ids, replace=True,
                                    size=(chromosomes.size, n_samples))

    # Tidy chromosome-swapped samples
    chromosomes_idx = pd.Index(chromosomes, name="Chromosome")
    samples_swapped_df = pd.DataFrame(samples_swapped_np, index=chromosomes_idx)
    samples_swapped_df = (samples_swapped_df
                          .melt(value_name="sample_id",
                                var_name="sample_id_sim",
                                ignore_index=False)
                          .reset_index())

    # Prepare data to match template
    data = (samples_swapped_df
            .assign(sample_id_sim=\
                simulate_sample_ids(samples_swapped_df["sample_id_sim"]))
            .merge(template_copynumber)
            .drop("sample_id", axis=1)
            .rename(columns=dict(sample_id_sim="sample_id")))
    
    return data

def simple_copynumber(n_samples, n_vars, *args, random_state=None, **kwargs):
    copynumber = simulate_matrix(n_samples, n_vars, "cnvr_", random_state)
    copynumber = copynumber.rename_axis("sample_id", axis=0).rename_axis("cnvr_id", axis=1)
    return copynumber

def simple_copynumber_segments(n_samples, *args, **kwargs):
    data_per_chrom = []
    for chrom, chrom_size in tqdm.tqdm(CHROMOSOME_SIZES.items()):
        chrom_data = simple_copynumber_segments_per_chrom(
            n_samples,
            genome_size = chrom_size,
            **kwargs
        )
        chrom_data = chrom_data.assign(Chromosome=chrom)
        data_per_chrom.append(chrom_data)
    data = pd.concat(data_per_chrom, axis=0)
    return data

def simple_copynumber_segments_per_chrom(
    n_samples,
    genome_size = 1e8,
    n_breakpoints_per_mb = 10,
    frac_breakpoints_remove = 0.6,
    diploid_val = 0,
    random_state=None,
):
    rng = np.random.default_rng(random_state)
    n_breakpoints = int(genome_size / 1e6 // n_breakpoints_per_mb)

    # Make breakpoints
    breakpoints = rng.integers(0, genome_size,
                               size=(n_breakpoints-1, n_samples))
    breakpoints = np.r_[ np.zeros((1, n_samples)), breakpoints ]
    breakpoint_remove = rng.choice([0, 1],
                                   p=[1-frac_breakpoints_remove,
                                      frac_breakpoints_remove],
                                   size=breakpoints.shape
                                   ).astype(bool)

    # Segment starts
    breakpoints_start = np.sort(breakpoints, axis=0)
    breakpoints_start[breakpoint_remove] = np.nan
    breakpoints_start[0, :] = 0. # Replace zeros that were removed
    breakpoints_start = np.sort(breakpoints_start, axis=0)

    # Segment ends
    breakpoints_end = np.roll(breakpoints_start, -1, axis=0) - 1
    breakpoints_end[-1, :] = genome_size
    breakpoints_end = np.sort(breakpoints_end, axis=0)

    # Segment values
    copynumber_vals = np.full(shape=breakpoints.shape,
                              fill_value=diploid_val, dtype=float)
    copynumber_vals[breakpoint_remove] = np.nan
    copynumber_vals = np.sort(copynumber_vals, axis=0)

    # Simulate SCNAs
    copynumber_vals = rng.normal(copynumber_vals, 1)
    
    copynumber_dict = []
    for i in range(n_breakpoints):
        for j in range(n_samples):
            copynumber_val = copynumber_vals[i, j]
            if pd.isna(copynumber_val):
                continue
            else:
                copynumber_dict.append(
                    dict(Start=breakpoints_start[i, j],
                        End=breakpoints_end[i, j],
                        sample_id=simulate_sample_ids([j])[0],
                        value=copynumber_val)
                )

    data = pd.DataFrame.from_dict(copynumber_dict)
    data["Start"] = data["Start"].astype(int)
    data["End"] = data["End"].astype(int)
    
    return data

def take_copynumber_subset(n_samples, *args, template_copynumber: pd.DataFrame = None, **kwargs):
    return take_subset(n_samples, template_data=template_copynumber, **kwargs)

copynumber_methods = {
    "simple_copynumber": simple_copynumber,
    "simple_copynumber_segments": simple_copynumber_segments,
    "take_copynumber_subset": take_copynumber_subset,
    "chromswap_augmented": chromswap_augmented
}