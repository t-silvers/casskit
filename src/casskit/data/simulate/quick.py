from dataclasses import dataclass
from typing import List, Tuple

import pandas as pd

import casskit.data.simulate.base as base


@dataclass
class QuickSim(base.SimulationMixin):    
    sampling_dist: str
    shape: Tuple
    dims: List
    stems: List = None

    @property
    def coord_labels(self):
        return [
            self.annotate(
                self.stems[__], self.shape[__]
            ) for __ in range(len(self.shape))
        ]

    @property
    def function_aliases(self):
        return dict(
            zeros=self.zeros,
            integers=lambda x: self.integers(size=x),
        )

    def check_inputs(self):
        if len(self.dims) != len(self.shape):
            raise ValueError("len(dims) != len(shape)")

        if self.stems is not None:
            if len(self.dims) != len(self.stems):
                raise ValueError("len(dims) != len(stems)")

        else:
            self.stems = self.dims

    def simulate(self) -> pd.DataFrame:
        data = self.function_aliases[self.sampling_dist](self.shape)

        return pd.DataFrame(data,
                            index=self.coord_labels[0],
                            columns=self.coord_labels[1])

    @classmethod
    def from_dict(cls, d) -> pd.DataFrame:
        return cls(**d).data

    def __post_init__(self):
        super().__init__()
        self.check_inputs()
        self.data = self.simulate()


quick_sim = QuickSim.from_dict
""""Convenience function for creating a quick simulation."""