# from dask_ml import compose as dask_compose
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.compose import ColumnTransformer, make_column_selector

__all__ = ["ColumnSelector", "LHSelector"]


class LHSelector(BaseEstimator, TransformerMixin):
    """Transformer selector for model LHS."""
    def __init__(self, col, samples):
        self.col = col
        self.samples = samples

    def fit(self, X, y=None):
        return self

    def fit_transform(self, X, y=None):
        """Override mixin fit_transform."""
        return self.transform(X)

    def transform(self, X, y=None):
        # TODO: Imputer
        return (X.filter(self.col)
                .reindex(self.samples)
                .fillna(0))

class ColumnSelector(ColumnTransformer):
    # For dask objects, we need to use the dask-ml version of ColumnTransformer
    # https://ml.dask.org/modules/generted/dask_ml.compose.ColumnTransformer.html#dask_ml.compose.ColumnTransformer
    # dask_ml.compose.ColumnTransformer(preserve_dataframe=False)
    def __init__(self, selector_name, **selector_kwargs):
        self.selector_name = selector_name
        super().__init__(
            transformers=[(self.selector_name, "passthrough",
                           make_column_selector(**selector_kwargs))],
            remainder="drop"
        )

# TODO:
# class ColumnSelectorDask(dask_compose.ColumnTransformer):
#     # For dask objects, we need to use the dask-ml version of ColumnTransformer
#     def __init__(self, selector_name, **selector_kwargs):
#         self.selector_name = selector_name
#         super().__init__(
#             transformers=[(self.selector_name, "passthrough",
#                            make_column_selector(**selector_kwargs))],
#             remainder="drop",
#             preserve_dataframe=False
#         )

class make_column_selector:
    """Create a callable to select columns to be used with
    :class:`ColumnTransformer`.

    :func:`make_column_selector` can select columns based on datatype or the
    columns name with a regex. When using multiple selection criteria, **all**
    criteria must match for a column to be selected.

    Parameters
    ----------
    pattern : str, default=None
        Name of columns containing this regex pattern will be included. If
        None, column selection will not be selected based on pattern.

    dtype_include : column dtype or list of column dtypes, default=None
        A selection of dtypes to include. For more details, see
        :meth:`pandas.DataFrame.select_dtypes`.

    dtype_exclude : column dtype or list of column dtypes, default=None
        A selection of dtypes to exclude. For more details, see
        :meth:`pandas.DataFrame.select_dtypes`.

    Returns
    -------
    selector : callable
        Callable for column selection to be used by a
        :class:`ColumnTransformer`.

    See Also
    --------
    ColumnTransformer : Class that allows combining the
        outputs of multiple transformer objects used on column subsets
        of the data into a single feature space.

    Examples
    --------
    >>> from sklearn.preprocessing import StandardScaler, OneHotEncoder
    >>> from sklearn.compose import make_column_transformer
    >>> from sklearn.compose import make_column_selector
    >>> import numpy as np
    >>> import pandas as pd  # doctest: +SKIP
    >>> X = pd.DataFrame({'city': ['London', 'London', 'Paris', 'Sallisaw'],
    ...                   'rating': [5, 3, 4, 5]})  # doctest: +SKIP
    >>> ct = make_column_transformer(
    ...       (StandardScaler(),
    ...        make_column_selector(dtype_include=np.number)),  # rating
    ...       (OneHotEncoder(),
    ...        make_column_selector(dtype_include=object)))  # city
    >>> ct.fit_transform(X)  # doctest: +SKIP
    array([[ 0.90453403,  1.        ,  0.        ,  0.        ],
           [-1.50755672,  1.        ,  0.        ,  0.        ],
           [-0.30151134,  0.        ,  1.        ,  0.        ],
           [ 0.90453403,  0.        ,  0.        ,  1.        ]])
    """

    def __init__(self, pattern=None, *, dtype_include=None, dtype_exclude=None):
        self.pattern = pattern
        self.dtype_include = dtype_include
        self.dtype_exclude = dtype_exclude

    def __call__(self, df):
        """Callable for column selection to be used by a
        :class:`ColumnTransformer`.
        
        Parameters
        ----------
        df : dataframe of shape (n_features, n_samples)
            DataFrame to select columns from.
        """
        if not hasattr(df, "iloc"):
            raise ValueError(
                "make_column_selector can only be applied to pandas dataframes"
            )
        df_row = df.iloc[:1]
        if self.dtype_include is not None or self.dtype_exclude is not None:
            df_row = df_row.select_dtypes(
                include=self.dtype_include, exclude=self.dtype_exclude
            )
        cols = df_row.columns
        if self.pattern is not None:
            cols = cols[cols.str.contains(self.pattern, regex=True)]
        return cols.tolist()