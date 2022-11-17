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

    ID: str
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
        return (f"{self.__class__.__name__}("
                f"name={self.ID}, "
                f"coords=({self.coords.chrom}, {self.coords.start_pos}, {self.coords.start_pos}), "
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
    
    k_cis: int = 2
    k_trans: int = 10
    k_var: int = 5

    var_ids: List = field(init=True, default=None)
    copynumber_ids: Dict = field(init=True, default=None)

    @property
    def cis_cn(self) -> List:
        return [Regulator(
            ID=self.egene_id, coords=self.egene_coords, cis=True, etype="cn"
        )]

    @property
    def var_eqtls(self) -> List:
        var_eqtls = self.rng.choice(self.var_ids, size=self.k_var, replace=False)
        return [
            Regulator(
                ID=__, cis=True, coords=self.egene_coords, etype="variant", efunc="linear"
            ) for __ in var_eqtls
        ]

    @property
    def cis_eqtls(self) -> List:
        cis_eqtls = []
        for _ in range(int(self.k_cis)):
            coords = self.randcncoords(chroms=[self.egene_coords.chrom])
            cis_eqtls.append(
                Regulator(ID="", cis=True, coords=coords, etype="variant", efunc="linear")
            )
            
        return cis_eqtls

    @property
    def trans_eqtls(self) -> List:
        trans_eqtls = []
        for _ in range(int(self.k_trans)):
            coords = self.randcncoords(chroms=list(self.copynumber_ids.keys()))
            trans_eqtls.append(
                Regulator(ID="", cis=False, coords=coords, etype="variant", efunc="linear")
            )
            
        return trans_eqtls

    def randcncoords(self, chroms, size=5E4) -> Coords:
        randchrom = self.rng.choice(chroms)
        randstart = self.rng.choice(self.copynumber_ids[randchrom])
        return Coords(randchrom, randstart, randstart+size)

    def __post_init__(self):
        super().__init__()
        if self.var_ids is None:
            self.var_ids = self.annotate("snp", self.k_var)

        if self.copynumber_ids is None:
            self.copynumber_ids = {k:v for k,v in zip(self.chroms, np.arange(0, 25E7, 1E6, dtype=int))}
        
        self.trans_chroms = list(set(self.chroms) - set([self.egene_coords.chrom]))
        self.regs = self.cis_cn + self.var_eqtls + self.cis_eqtls + self.trans_eqtls

    def df(self):
        return (
            pd.DataFrame.from_records([__.coords for __ in self.regs], columns=Coords._fields)
            .assign(is_cis=[__.cis for __ in self.regs],
                    feature_id=[__.ID for __ in self.regs],
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
        k_cis=2,
        k_trans=10,
        k_var=5,
        var_ids=None,
        copynumber_ids=None,
        return_df=False,
    ):
        egene_coords = Coords(*egene_coords)
        grn = cls(egene_id, egene_name, egene_coords, chroms, k_cis,
                  k_trans, k_var, var_ids, copynumber_ids)
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