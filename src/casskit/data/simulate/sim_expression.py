from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd
import scipy as sp

import casskit.data.simulate.base as base
import casskit.data.simulate.sim_grn as grn


class SimExpression(base.SimulationMixin):

    _data = None

    def __init__(
        self,
        gene_id,
        N: int = 100,
        p: int = 500,
        regulators: List[grn.Regulator] = None,
        variants: pd.DataFrame = None,
        copynumber: pd.DataFrame = None,
        exp_method="gauss",
        noise_sd: int=1,
        **kwargs
    ) -> None:
        super().__init__(N, p)
        self.gene_id = gene_id
        self.regulators = regulators
        self.variants = variants
        self.copynumber = copynumber
        self.exp_method = exp_method
        self.noise_sd = noise_sd
        self.kwargs = kwargs
        
        # Methods
        self.methods = {"gauss": self.gauss,
                        "poisson": self.poisson,
                        "linreg": self.linreg,
                        "poireg": self.poireg}

    @property
    def data(self):
        µ = self.sim_loc()
        y = µ + self.noise(self.noise_sd, self.N)
        return pd.DataFrame(y, index=self.annotate("TCGA", size=self.N),
                            columns=[f"expression_{self.gene_id}"])

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

    def sim_loc(self):
        µ = 0
        for eqtl in self.regulators:
            if eqtl.etype == "copynumber":
                # Note that this only honors the start position
                x = (self.copynumber
                     .query("""
                            Chrom == @eqtl.coords.chrom & \
                            Start <= @eqtl.coords.start_pos & \
                            End >= @eqtl.coords.start_pos
                            """)
                     .groupby("sample")
                     ["value"].mean().values)
                
            elif eqtl.etype == "variant":
                x = self.variants[eqtl.ID].values
            
            else:
                raise ValueError("Unknown eQTL type.")
            
            µ += eqtl.expression_contribution(x)

        return µ

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
    def simulate(
        cls,
        gene_id,
        N: int = 100,
        p: int = 500,
        regulators: List[grn.Regulator] = None,
        variants: pd.DataFrame = None,
        copynumber: pd.DataFrame = None,
        exp_method="gauss",
        **kwargs
    ) -> pd.DataFrame:
        return cls(gene_id, N, p, regulators, variants, copynumber, exp_method, **kwargs).data

simulate_expression = SimExpression.simulate
