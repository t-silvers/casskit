from collections.abc import Iterable

import pandas as pd


def simulate_ids(num, prefix=None, as_idx=False, idx_name=None):
    if isinstance(num, int):
        num = range(num)
    elif isinstance(num, Iterable):
        pass
    if prefix is None:
        sim_ids = num
    else:
        sim_ids = [f"{prefix}{i:04}" for i in num]
    if as_idx is False:
        return sim_ids
    else:
        if idx_name is None:
            idx_name = f"{prefix}id"
        return pd.Index(sim_ids, name=idx_name)

def simulate_sample_ids(num, **kwargs):
    return simulate_ids(num, prefix="sample_", **kwargs)