from collections import namedtuple
import string
from typing import Union
import warnings

import numpy as np
import pandas as pd
import xarray as xr


__all__ = ["MultiOmicData"]


class MultiOmicData:
    """Lightweight multi-omic class.
    
    Goals to add data type flexibility, advanced caching,
    advanced slicing, pyranges joins, etcc to openomic.MultiOmicData.
    """
    pass