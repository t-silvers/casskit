from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
import scipy as sp

import casskit.data.simulate.base as base


class SimExpression(base.SimulationMixin):

    _data = None

    def __init__(
        self,
        N: int = 100,
        p: int = 500,
        exp_method="gauss",
        **kwargs
    ) -> None:
        super().__init__(N, p)
        self.exp_method = exp_method
        self.kwargs = kwargs
        
        # Methods
        self.methods = {"gauss": self.gauss,
                        "poisson": self.poisson,
                        "linreg": self.linreg,
                        "poireg": self.poireg}

    @property
    def data(self):
        return self.make_data()

    def make_data(self):
        self.methods[self.cn_method](self.kwargs)

    def gauss(self):
        # ill-conditioned
        return self.rng.standard_normal((self.N, self.p))

    def poisson(self):
        lambd = self.kwargs.get("lambd", 1)
        return self.rng.poisson(lambd, size=(self.N, self.p))

    def linreg(self):
        noise_sd = self.kwargs.get("noise_sd", 1)
        intercept = self.kwargs.get("intercept", 1)
        design = self.kwargs.get("design", None)
        if design is None:
            raise ValueError("Must provide a design matrix.")
        
        betas = self.kwargs.get("betas", 1)
        
        assert len(betas) == design.shape[0], "Number of betas must equal number of features."
        
        ε = self.noise(noise_sd, self.N)
        β0 = self.intercept(intercept)
        µ = β0 + np.dot(design, betas) + ε
        
        return µ

    def poireg(self):
        µ = self.linreg()
        return np.exp(µ)