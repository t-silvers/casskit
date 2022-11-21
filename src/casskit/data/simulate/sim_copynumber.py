# Author: Thomas R. Silvers <thomas.silvers.1@gmail.com>
# License: MIT

from typing import List

import numpy as np
import pandas as pd
import ruptures as rpt

import casskit.data.simulate.base as base


def bamgineer():
    # https://github.com/pughlab/bamgineer
    pass

def xome_blender():
    # https://github.com/rsemeraro/XomeBlender
    pass

def cic():
    # Change in covariance
    # https://github.com/mlondschien/changeforest/tree/main/changeforest-py
    # https://github.com/mlondschien/changeforest-simulations
    pass


def nhpp_sim_constwindow():
    """Simulate a Non-Homogeneous PP with constant window spike
    
    Implements seqCBS::nhppSimConstWindowGen() from the
    [seqCBS R package](https://cran.r-project.org/web/packages/seqCBS/index.html)
    [paper](https://arxiv.org/abs/1206.6627)
    """
    pass


def nhpp_sim_const_window_gen(
    controlRates,
    filename,
    chromosomeN,
    nSpike=25,
    cptLen: List = [3,5,8,12,20,30,50,75,100],
    nPair=2,
    nRepeat=10,
    minGain=1.5,
    maxGain=4,
    minLoss=0.01,
    maxLoss=0.5,
    pGain=0.6,
):
    """Simulate a Non-Homogeneous PP with constant window spike
    
    Implements seqCBS::nhppSimConstWindowGen() from the
    [seqCBS R package](https://cran.r-project.org/web/packages/seqCBS/index.html)
    [paper](https://arxiv.org/abs/1206.6627)
    """
    pass


def pw_constant(n_samples: int = 500, dim: int = 3, n_bkps: int = 6, sigma: float = 1.0) -> np.ndarray:
    """Piecewise constant signal.
    
    Thin wrapper around ruptures.datasets.pw_constant.

    Args:
        n_samples (int): number of samples
        dim (int): dimension
        n_bkps (int): number of breakpoints
        sigma (float): standard deviation of the noise

    Returns:
        np.ndarray: signal
        np.ndarray: breakpoints
    """
    bkps = np.sort(np.random.randint(0, n_samples, size=n_bkps))
    signal, bkps = rpt.pw_constant(n_samples, dim, n_bkps, noise_std=sigma)

    return signal, bkps


def chrom_swap_augmented(
    copynumber: pd.DataFrame,
    n_samples: int = 500,
    rng = np.random.default_rng(),
):
    samples = copynumber["sample"].unique()
    groups = copynumber["Chrom"].unique()

    return (
        pd.DataFrame(rng.choice(samples,
                                replace=True,
                                size=[len(groups), n_samples]),
                    index=groups)
        .rename_axis("Chrom")
        .melt(value_name="sample",
              var_name="sample_sim",
              ignore_index=False)
        .reset_index()
        .merge(copynumber)
        .drop("sample", axis=1)
        .rename(columns={"sample_sim": "sample"})
    )

class SimCopynumber(base.SimulationMixin):
    """Simulated copynumber data.

    Simulated copy number by various method.

    Typical usage example:

    simulated_copynumber = SimCopynumber.simulate(cn_method="gauss")
    """

    _data = None

    def __init__(
        self,
        N: int = 100,
        p: int = 500,
        cn_method="gauss",
        **kwargs
    ) -> None:
        super().__init__(N, p)
        self.cn_method = cn_method
        self.kwargs = kwargs
        
        # Methods
        self.methods = {"gauss": self.gauss, "markov": self.markov,
                        "swap_augment": self.swap_augment,
                        "original": self.original}

    @property
    def data(self):
        return self.make_data()

    def make_data(self):
        sim_data = self.methods[self.cn_method]()
        sim_data["sample"] = self.annotate("TCGA-00", range_=sim_data.pop("sample"))
        return sim_data

    def gauss(self):
        # ill-conditioned
        return self.rng.standard_normal((self.N, self.p))

    def original(self):
        return self.kwargs.get("copynumber_template", None) 

    def markov(self):
        G = self.kwargs.get("groups", 1)
        X = np.zeros([self.N, self.p])
        if G > 1:
            X_ = []
            for X_g in np.array_split(X, G, axis=1):
                X_.append(self._x_mchain_w_restarts(X_g))
            return np.concatenate(X_, axis=1)
        
        else:
            return self._x_mchain_w_restarts(X)

    def _x_mchain_w_restarts(
        self,
        X,
        transition_probs: List[float] = [0.05, 0.9, 0.05],
        restart_prob: float = 0.1
    ):
        """Simulate SCNA data as a Markov chain with restarts.
        
        Note: There are no boundary conditions, so values can, for instance,
        dip to zero and lower while being subject to the same transition
        probabilities.
        """
        np.testing.assert_almost_equal(sum(transition_probs), 1)
        
        X[:,0] = 2
        for i in range(1, X.shape[1]):
            X[:,i] = X[:,i-1] + self.rng.choice([-1, 0, 1], p=transition_probs, size=X.shape[0])
            if self.rng.binomial(1, restart_prob):
                X[:,i] = 2
        
        return X

    def swap_augment(self):
        """Data augmentation by swapping chromomes.
        
        Args:
            copynumber (pd.DataFrame): segment copynumber data
        
        Returns:
            augmented_ix: pd.DataFrame
                Index of augmented samples.
                | Chrom || sample || Start || End || value |
                where sample is the simulated sample
                
        """
        copynumber_template = self.kwargs.get("copynumber_template", None) 
               
        return chrom_swap_augmented(
            copynumber=copynumber_template,
            n_samples=self.N,
            rng=self.rng,
        )

    def obs_X(self, X, X_sd=0.1):
        X += self.rng.normal(0, X_sd, size=X.shape)
        X[X < 0] = 0
        return X

    @classmethod
    def simulate(
        cls,
        N: int = 100,
        p: int = 500,
        cn_method="gauss",
        **kwargs
    ) -> None:
        sim = cls(N, p, cn_method, **kwargs)
        return sim.data
    

simulate_copynumber = SimCopynumber.simulate
"""Simulate copynumber data.

Args:
    n_samples (int): number of samples
    dim (int): dimension
    n_bkps (int): number of breakpoints
    sigma (float): standard deviation of the noise

Returns:
    np.ndarray: signal
    np.ndarray: breakpoints

Examples:
    By augmenting data:

    >>> import casskit as ck
    >>> copynumber_raw = ck.io.get_tcga("cnv", cancer)
    >>> ck.data.simulate_copynumber(cn_method="swap_augment", copynumber_template=copynumber_raw)

"""