from collections import namedtuple
from dataclasses import dataclass, field
from typing import Any, Dict, Callable, List, Tuple

import numpy as np
import pandas as pd

import casskit.data.simulate.base as base


Coords = namedtuple('Coords', ['chrom', 'start_pos', 'end_pos'])

@dataclass
class Regulator(base.SimulationMixin):
    """An eGene regulator."""

    efunc_names = ['log', 'log2', 'tanh']

    coords: Coords
    cis: bool = True
    etype: str = "variant"
    efunc: str = "linear"

    @property
    def functional_relationship(self) -> Callable:
        """The functional relationship of the regulator."""
        if self.efunc == "linear":
            return lambda x: x
        return getattr(np, self.efunc)

    @property
    def beta(self) -> float:
        return self.rng.normal(loc={True: 2, False: 0}[self.cis], scale=1)

    def expression_contribution(self, x):
        return self.functional_relationship(x) * self.beta

    def __repr__(self):
        return (f"{self.__class__.__name__} "
                f"(coords=({self.coords.chrom}, {self.coords.start_pos}, {self.coords.start_pos}), "
                f"is_cis={self.cis}, "
                f"eqtl_type={self.etype}, "
                f"eqtl_func={self.efunc})")

    def __str__(self):
        return self.__repr__()
    
    def __post_init__(self):
        super().__init__()


@dataclass
class SimGRN(base.SimulationMixin):

    egene_id: str
    egene_name: str
    egene_coords: Coords

    chroms: List[str] = None
    trans_chroms: List[str] = field(init=False, default=None)
    
    k_cis: int = 5
    k_trans: int = 10

    @property
    def cis_cn(self) -> List:
        return [Regulator(coords=self.egene_coords, cis=True, etype="cn")]

    @property
    def cis_eqtls(self) -> List:
        cis_eqtls = []
        for _ in range(int(self.k_cis)):
            coords = self.randcoords(chroms=[self.egene_coords.chrom],
                                     starts=np.arange(0, 25E7, 1E6, dtype=int))
            cis_eqtls.append(
                Regulator(cis=True, coords=coords, etype="variant", efunc="linear")
            )
            
        return cis_eqtls

    @property
    def trans_eqtls(self) -> List:
        trans_eqtls = []
        for _ in range(int(self.k_trans)):
            coords = self.randcoords(chroms=self.trans_chroms,
                                     starts=np.arange(0, 25E7, 1E6, dtype=int))
            trans_eqtls.append(
                Regulator(cis=False, coords=coords, etype="variant", efunc="linear")
            )
            
        return trans_eqtls

    def randcoords(self, chroms, starts, size=1) -> Coords:
        randchrom = self.rng.choice(chroms)
        randstart = self.rng.choice(starts)
        return Coords(randchrom, randstart, randstart+size)

    def __post_init__(self):
        super().__init__()
        self.trans_chroms = list(set(self.chroms) - set([self.egene_coords.chrom]))
        self.regs = self.cis_cn + self.cis_eqtls + self.trans_eqtls

    def df(self):
        return (
            pd.DataFrame.from_records([__.coords for __ in self.regs], columns=Coords._fields)
            .assign(is_cis=[__.cis for __ in self.regs],
                    eqtl_type=[__.etype for __ in self.regs],
                    eqtl_func=[__.efunc for __ in self.regs],
                    beta=[__.beta for __ in self.regs])
        )

    @classmethod
    def simulate(
        cls,
        egene_id,
        egene_name,
        egene_coords,
        chroms,
        k_cis=5,
        k_trans=10,
        return_df=False,
    ):
        egene_coords = Coords(*egene_coords)
        grn = cls(egene_id, egene_name, egene_coords, chroms, k_cis, k_trans)
        if return_df is True:
            return grn.regs, grn.df()
        return grn.regs

simulate_grn = SimGRN.simulate
"""Simulate gene regulatory network data.

Args:
    egene info

Returns:
    pd.DataFrame

Examples:
    By augmenting data:

    >>> import casskit as ck
    >>> copynumber_raw = ck.io.get_tcga("cnv", cancer)
    >>> ck.data.simulate_grn(
        egene_id='ENSG000001', egene_name='EGFR', egene_coords=Coords(chrom='chr7', start_pos=55000000, end_pos=55000001),
        chroms=['chr7', 'chr8', 'chr9'], trans_chroms=['chr9', 'chr8'], k_cis=5, k_trans=10
    )

"""