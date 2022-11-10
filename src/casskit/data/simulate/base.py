from typing import List

import numpy as np


def get_seed():
    return np.random.SeedSequence().entropy

class FunctionMixin:

    def intercept(self, a):
        return a
    
    def noise(self, sd, size):
        return self.rng.normal(0, sd, size=size)

    def linear(self, X, β):
        return np.dot(X, β)

    def poisson_reg(self, β0, X, β, σ):
        µ = self.intercept(β0) + self.linear(X, β) + self.noise(σ, X.shape[0])
        return np.exp(µ)

    def negbin_reg(self):
        pass


class SimulationMixin(FunctionMixin):
    
    _X = None
    _β = None
    _y = None
    
    def __init__(
        self,
        N=100,
        p=500,
        seed=None
    ) -> None:
        
        self.N = N
        self.p = p
        self.seed = seed if seed is not None else get_seed()
        self.rng = np.random.default_rng(self.seed)
    
    @staticmethod
    def annotate(stem, size) -> List[str]:
        return [f"{stem}_{i:04}" for i in range(size)]