import pandas as pd

import casskit.data.simulate.base as base


class SimVariants(base.SimulationMixin):

    _data = None

    def __init__(
        self,
        N: int = 100,
        p: int = 500,
        maf_method="binary",
        **kwargs
    ) -> None:
        super().__init__(N, p)
        self.maf_method = maf_method
        self.kwargs = kwargs
        
        # Methods
        self.methods = {"binary": self.binary}

    @property
    def data(self):
        self.annotate("snp", size=self.p)
        return pd.DataFrame(
            self.make_data(),
            columns=self.annotate("snp", size=self.p),
            index=self.annotate("TCGA-", size=self.N),
        )

    def make_data(self):
        return self.methods[self.maf_method](**self.kwargs)

    def binary(self, **kwargs):
        mut_freq = kwargs.get("mut_freq", 0.01)
        return self.rng.choice([0, 1],
                               size=(self.N, self.p),
                               p=[1 - mut_freq, mut_freq])

    @classmethod
    def simulate(cls, N=100, p=500, maf_method="binary", **kwargs):
        return cls(N, p, maf_method, **kwargs).data


simulate_variants = SimVariants.simulate
"""Simulate variants data.

Args:
    n_samples (int): number of samples

Returns:
    pd.DataFrame: simulated variants data

Examples:
    >>> import casskit as ck

"""