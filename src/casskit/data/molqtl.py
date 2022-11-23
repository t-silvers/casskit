# Class for dealing with molecular QTL data

from abc import ABC, abstractmethod

import pandas as pd

from casskit.typing import DATAFRAME


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
    

class cneQTL(MolQTL):
    def __init__(
        self,
        mol: str = "rna",
        qtl: str = "cn",
    ):
        super().__init__(mol, qtl)
    
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
    
