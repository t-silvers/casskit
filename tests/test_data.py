from pathlib import Path
import pytest

import pandas as pd

import casskit.data


def test_molqtl():
    
    # Test to check that casskit.data.molqtl checks column names
    with pytest.raises(ValueError):
        casskit.data.cneQTL().load_data(pd.DataFrame([1, 2, 3]))
    
    # Test to check that casskit.data.molqtl checks data type
    with pytest.raises(TypeError):
        casskit.data.cneQTL().load_data(pd.Series([1, 2, 3]))
    
    # Test to check that casskit.data.molqtl checks data type
    with pytest.raises(TypeError):
        casskit.data.cneQTL().load_data([1, 2, 3])
    
    # works with asset data
    data = read_eqtl_data()
    casskit.data.cneQTL().load_data(data)
    
    # check that works with random data with correct columns
    casskit.data.cneQTL().load_data(
        pd.DataFrame({
            "target": [1, 2, 3],
            "qtl": [1, 2, 3],
            "beta": [1, 2, 3],
            "proximity": [1, 2, 3],
        })
    )
    