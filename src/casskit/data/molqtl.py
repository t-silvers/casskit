# Class for dealing with molecular QTL data

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd

from casskit.typing import DATAFRAME


@dataclass
class Locus(ABC):
    chrom: str
    chrom_start: int
    chrom_end: int = None
    ID: str = None
    size: int = field(init=False)

    def __post_init__(self):
        if not self.chrom.startswith("chr"):
            self.chrom = f"chr{self.chrom}"
            
        if self.chrom_end is None:
            self.chrom_end = self.chrom_start

        # check that the start and end are ints
        self.chrom_start = int(self.chrom_start)
        self.chrom_end = int(self.chrom_end)

        # check that the start is less than the end
        if self.chrom_start > self.chrom_end:
            raise ValueError(
                f"chrom_start {self.chrom_start} is greater than chrom_end {self.chrom_end}"
            )
        
        self.size = max(self.chrom_end - self.chrom_start, 1)

@dataclass
class cneQTL(Locus): ...

class MolQTL(ABC):
    
    COLUMNS = [
        "target",
        "qtl",
        "beta",
        "proximity"
    ]
    
    def __init__(
        self,
        mol: str,
        qtl: str,
    ):
        self.mol = mol
        self.qtl = qtl

    def predict_effect(self):
        pass
    
    @abstractmethod
    def load_data(self):
        pass

    @abstractmethod
    def validate_data(self):
        pass

    @abstractmethod
    def is_cis(self):
        pass
    

class MolQTLDataFrame(MolQTL):
    pass

class MolQTLTarget:
    pass

class cneQTL(MolQTL):
    def __init__(
        self,
        mol: str = "rna",
        qtl: str = "cn",
    ):
        super().__init__(mol, qtl)
        
    def add(self, other):
        pass
    
    def load_data(self, data: DATAFRAME = None):
        self.data = self.validate_data(data)
        
        return self

    def validate_data(self):
        if not isinstance(self.data, pd.DataFrame):
            raise TypeError("Data must be a pandas DataFrame")

        if not self.data.columns.isin(self.COLUMNS).all():
            raise ValueError("Data must have columns: " + ", ".join(self.COLUMNS))

    def is_cis(self):
        # If on same chromosome, then cis
        pass


@dataclass(frozen=True)
class cneQTL:
    target: str
    qtl: str
    