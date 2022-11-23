import logging
from pathlib import Path
import pytest
import re

import pandas as pd


TEST_ASSETS = Path(__file__).parent / "assets"


def setup_logging(level):
    if level is None:
        level = logging.WARN
    
    elif re.match(r"\d+", level):
        level = int(level)
    
    logging.basicConfig(level=level)

@pytest.fixture
def data_dir(tmpdir_factory):
    """Setup base cache directory for a test"""
    p = tmpdir_factory.mktemp(".cache")
    
    return p

@pytest.fixture
def read_eqtl_data():
    return pd.read_csv(TEST_ASSETS / "eqtl_data.csv")