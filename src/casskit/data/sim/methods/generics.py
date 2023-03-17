import warnings

import numpy as np
import pandas as pd

from ..utils import simulate_ids, simulate_sample_ids


def simulate_matrix(n_samples, n_vars, var_prefix, random_state=None, rv="integers"):
    rng = np.random.default_rng(random_state)
    if rv == "integers":
        values = rng.integers(0, 1000, size=(n_samples, n_vars))
    elif rv == "normal":
        values = rng.normal(0, 1, size=(n_samples, n_vars))
    elif rv == "uniform":
        values = rng.uniform(0, 1, size=(n_samples, n_vars))
    df = pd.DataFrame(values,
                      index=simulate_sample_ids(n_samples),
                      columns=simulate_ids(n_vars, var_prefix))
    return df

def take_subset(
    n_samples,
    *args,
    subset_sample_ids=None, # May want to refactor elsewhere, otherwise n_samples should be kwarg
    template_data: pd.DataFrame = None,
    random_state=None,
    **kwargs,
):
    if subset_sample_ids is None:
        rng = np.random.default_rng(random_state)
        template_n_samples = template_data["sample_id"].nunique()
        
        # Check if there are enough samples in the template
        if template_n_samples < n_samples:
            warnings.warn(f"Requested {n_samples} samples, "
                        f"but only {template_n_samples} are available.")
            subset_sample_ids = template_data["sample_id"].unique()

        else:
            all_sample_ids = template_data["sample_id"].unique()
            subset_sample_ids = rng.choice(all_sample_ids, n_samples, False)
            
    template_data_subset = template_data.query("sample_id in @subset_sample_ids")
    sample_codes, __ = pd.factorize(template_data_subset["sample_id"])
    template_data_subset["sample_id"] = simulate_sample_ids(sample_codes)

    return template_data_subset