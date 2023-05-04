from typing import Callable

import sys
sys.path.append("/Users/thomassilvers/GitHub/")
from casskit_data import omic

import pandas as pd

from .base import OmicsBuilderMixin
from ..methods import simulation_methods
from ..utils import simulate_sample_ids, simulate_ids


class OmicsBuilder(OmicsBuilderMixin):
    def __init__(self, omic_type, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.omic_type = omic_type
        self.simulation_methods = simulation_methods[self.omic_type]

    @property
    def sim_data(self):
        return getattr(omic, self.omic_type)(data=self._data)

    @sim_data.setter
    def sim_data(self, data: pd.DataFrame):
        self._data = data

    @OmicsBuilderMixin.simulation_methods.setter
    def simulation_methods(self, methods):
        self._simulation_methods.update(methods)

    def simulate(self, method, *args, **kwargs):
        if isinstance(method, Callable):
            simulation_method = method
        else:
            simulation_method = self.simulation_methods.get(method, None)
            if simulation_method is None:
                raise ValueError(f"Method {method} not recognized")
        
        return simulation_method(
            self.n_samples,
            self.n_vars,
            *args,
            **kwargs
        )
    
class CopyNumberVariationBuilder(OmicsBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__("CopyNumberVariation", *args, **kwargs)

    def simulate_data(self, method, parse_output=True, *args, **kwargs):
        sim_data = self.simulate(method, *args, **kwargs)
        if parse_output is True:
            sample_idx = pd.Index(simulate_sample_ids(self.n_samples), name="sample_id")
            cnvr_idx = pd.Index(simulate_ids(self.n_vars, "cnvr_"), name="cnvr_id")
            self.sim_data = pd.DataFrame(sim_data, columns=cnvr_idx, index=sample_idx)
        else:
            self.sim_data = sim_data

class MessengerRNABuilder(OmicsBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__("MessengerRNA", *args, **kwargs)

    def simulate_data(self, method, *args, **kwargs):
        sim_data = self.simulate(method, *args, **kwargs)
        self.sim_data = sim_data

class SomaticMutationBuilder(OmicsBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__("SomaticMutation", *args, **kwargs)

    def simulate_data(self, method, *args, **kwargs):
        sim_data = self.simulate(method, *args, **kwargs)
        self.sim_data = sim_data

class MethylationBuilder(OmicsBuilder):
    def __init__(self, *args, **kwargs):
        super().__init__("Methylation", *args, **kwargs)

    def simulate_data(self, method, *args, **kwargs):
        sim_data = self.simulate(method, *args, **kwargs)
        self.sim_data = sim_data