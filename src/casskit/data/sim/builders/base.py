from abc import ABC, abstractmethod

import pandas as pd

from ..methods.generics import take_subset


class OmicSimulator:
    def __init__(self, n_samples, n_vars):
        self.n_samples = n_samples
        self.n_vars = n_vars
        self.builder = None
    
    def simulate_omics(self, builder, method, **kwargs):
        # Set attrs on simulator
        self.method = method
        self.method_kwargs = kwargs
        
        # Set attrs on builder
        self.builder = builder
        self.builder.n_samples = self.n_samples
        self.builder.n_vars = self.n_vars
        
        # Run builder
        self.builder.simulate_data(method, **kwargs)
    
    @property
    def simulated_omics(self):
        return self.builder.data

class OmicsBuilderMixin(ABC):
    _data = pd.DataFrame()
    _simulation_methods = {
        "take_subset": take_subset,
    }
    def __init__(
        self,
        *args,
        n_samples: int = None,
        n_vars: int = None,
        **kwargs,
    ):
        self.n_samples = n_samples
        self.n_vars = n_vars

    @property
    def simulation_methods(self):
        return self._simulation_methods

    @simulation_methods.setter
    @abstractmethod
    def simulation_methods(self):
        pass

    def simulate_data(self, method="simple_expression", *args, **kwargs):
        simulation_method = self.simulation_methods.get(method, None)
        if simulation_method is None:
            raise ValueError(f"Method {method} not recognized")
        self.sim_data = simulation_method(self.n_samples, self.n_vars, **kwargs)