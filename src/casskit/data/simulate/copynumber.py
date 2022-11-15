from typing import List

import numpy as np
import ruptures as rpt

import casskit.data.simulate.base as base


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

class SimCopynumber(base.SimulationMixin):

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
        self.methods = {"gauss": self.gauss,
                        "markov": self.markov}

    @property
    def data(self):
        return self.make_data()        

    def make_data(self):
        self.methods[self.cn_method]()

    def gauss(self):
        # ill-conditioned
        return self.rng.standard_normal((self.N, self.p))

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

    def obs_X(self, X, X_sd=0.1):
        X += self.rng.normal(0, X_sd, size=X.shape)
        X[X < 0] = 0
        return X