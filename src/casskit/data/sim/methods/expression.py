import numpy as np
import pandas as pd

from .generics import simulate_matrix


def from_grn(
    n_samples: int, # TODO: Unused
    *args,
    grn: pd.DataFrame = None,
    design: pd.DataFrame = None,
    noise_sigma: float = 1,
    rng = np.random.default_rng(),
    **kwargs,
):
    # TODO: Validators
    # TODO: Make more robust
    expression_det = design.mul(grn["beta"], axis="index").groupby("gene_id").sum()
    expression = expression_det \
        + rng.normal(0, noise_sigma, size=expression_det.shape)
    expression = expression.T
    return expression

def simple_expression(n_samples, n_vars, random_state=None, **kwargs):
    expression = simulate_matrix(n_samples, n_vars, "gene_", random_state)
    expression = expression.rename_axis("sample_id", axis=0).rename_axis("gene_id", axis=1)
    return expression

expression_methods = {
    "simple_expression": simple_expression,
    "from_grn": from_grn,
}