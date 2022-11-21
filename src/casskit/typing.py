from typing import Dict, List, Tuple, Union

import dask.dataframe as dd
import dask.distributed
import pandas as pd

DATAFRAME = Union[dd.DataFrame, pd.DataFrame, dask.distributed.Future]