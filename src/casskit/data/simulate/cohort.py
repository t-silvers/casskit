from pathlib import Path

import pandas as pd

import casskit.data.simulate.base as base
from casskit.data.simulate.copynumber import SimCopynumber
from casskit.data.simulate.expression import SimExpression
from casskit.data.simulate.variants import SimVariants


class SimEQTLCohort(base.SimulationMixin):
    def __init__(self, N=100, seed=None, **kwargs) -> None:
        super().__init__(N, seed)
        self.kwargs = kwargs

    @property
    def betas(self):
        # Simulate betas
        pass
    
    @property
    def data(self):
        return self.make_data()

    @property
    def design(self):
        return self.make_design()
    
    @property
    def copynumber(self):
        cn_p = self.kwargs.get("cn_p", 500)
        cn_method = self.kwargs.get("cn_method", "markov")
        
        return SimCopynumber(self.N, cn_p, cn_method, **self.kwargs).data

    @property
    def expression(self):
        ex_p = self.kwargs.get("ex_p", 500)
        exp_method = self.kwargs.get("exp_method", "linreg")
        exp_params = self.kwargs.update(dict(betas=self.betas, design=self.design))
        
        return SimExpression(self.N, ex_p, exp_method, **exp_params).data

    @property
    def variants(self):
        maf_p = self.kwargs.get("maf_p", 500)
        return SimVariants(self.N, maf_p, **self.kwargs).data

    def make_data(self):
        # Combine copynumber and variants and expression and betas
        # Annotate with labels
        pass

    def make_design(self):
        # Combine copynumber and variants
        pass