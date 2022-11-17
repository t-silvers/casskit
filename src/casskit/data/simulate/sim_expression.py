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

    @property
    def latent_vars(self) -> List:
        # https://cran.r-project.org/web/packages/simstandard/vignettes/simstandard_tutorial.html
        # https://bookdown.org/marklhc/notes/simulation-example-on-structural-equation-modeling-sem.html
        
        alphas = self.rng.normal([0, 0.5])
        phi = np.array([[1, 0.1], [0.1, 0.2]])
        lambd = np.array([[1, 1, 1, 1], [0, 1, 2, 3]])
        theta = np.diag([0.5, 0.5, 0.5, 0.5])
        eta = self.rng.multivariate_normal(mean=alphas, cov=phi, size=self.N)
        e = self.rng.multivariate_normal(mean=[0, 0, 0, 0], cov=theta, size=self.N)
        y = np.dot(eta, lambd) + e
        
        return pd.DataFrame(y)

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

    @classmethod
    def linmod(cls, X, link_func=lambda x: x):
        pass


variant_sim
copynumber_sim

grn_sim[0].expression_contribution(-4)



simulate_expression = SimExpression.linmod