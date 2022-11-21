from typing import List

import numpy as np


# https://github.com/yaglm/yaglm/blob/bf8169aa73ae1ecea8f8bb0f6d3c093886dc8507/yaglm/toy_data.py


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
    def annotate(stem, size: int = None, range_ : List = []) -> List[str]:
        if size is not None:
            # return map(lambda x: f"{stem}_{x:04}", range(size))
            return [f"{stem}-{i:04}" for i in range(size)]
        elif len(range_) > 0:
            return [f"{stem}-{i:04}" for i in range_]
        else:
            raise ValueError("Must provide either size or range_.")
    