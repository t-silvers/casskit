from typing import Dict, List, Tuple, Union

import dask.dataframe as dd
from dask.distibuted import Future
import pandas.DataFrame as pd


DATAFRAME = Union[dd.DataFrame, pd.DataFrame, Future]