import imp
from typing import Tuple
import pandas as pd


from casskit.typing import DATAFRAME


def feature_selector(
    model_frame: DATAFRAME,
    ixs: tuple(int, range),
    ret_df: bool = True,
    fillna: bool = True,
):
    response_ix, design_ix = ixs
    
    X = model_frame.iloc[:, design_ix]
    if fillna:
        X = X.fillna(X.mean())

    y = model_frame.iloc[:, response_ix]

    if ret_df:
        return X, y
    return X.values, y.values
