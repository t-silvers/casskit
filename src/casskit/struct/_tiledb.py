from ast import Expression
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from pathlib import Path
import shutil
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import tiledb
from tiledb.libtiledb import TileDBError

from casskit import __version__
from ..config import CACHE_DIR # TEMP

cfg = tiledb.Ctx().config()
cfg.update(
  {
    'py.init_buffer_bytes': 1024**2 * 50
  }
)
tiledb.default_ctx(cfg)


class DataDB:
    """Turns a pandas dataframe into a TileDB array"""
    def __init__(
        self,
        data: pd.DataFrame,
        feature: str,
        tiledb_grp: Dict[str, str],
        cache_dir: Path = CACHE_DIR,
        overwrite: bool = False
    ) -> None:
        self.data = data
        self.feature = feature
        self.tiledb_grp = tiledb_grp
        self.cache_dir = cache_dir
        self.overwrite = overwrite
        
        # Make db
        print("Making TileDB database...")
        self._make_db()

    def _make_db(self) -> None:
        self.uris = self._make_uris(feature=self.feature)
        if (not self.uris["data"].exists() or self.overwrite is True):
            self._make_data_array()
            self._make_label_arrays()
            self._ingest()

    @staticmethod
    def _meta() -> Dict[str, str]:
        return f"casskit_version {__version__}"
    
    @staticmethod
    def _uri(cache_dir: Path, groups: Dict[str, str], dbtype: str) -> Path:
        uri = Path(cache_dir)
        for k,v in groups.items():
            uri /= f"{k}~{v}"
            try:
                tiledb.group_create(f"{k}~{v}")
            except TileDBError as err:
                print("TileDB exception: ", err)
            
        uri /= f"{dbtype}.tiledb"
        uri.parent.mkdir(parents=True, exist_ok=True)
        
        return uri

    def _make_uris(self, **kwargs) -> Dict[str, Path]:
        """Make uris for data and labels"""
        uris = {}
        for __ in ["data", "sample", "feature"]:
            uris[__] = self._uri(self.cache_dir,
                                 self.tiledb_grp,
                                 dbtype=__)
            
            if (uris[__].exists() and self.overwrite is True):
                shutil.rmtree(uris[__])
        
        return uris

    def _make_data_array(self) -> None:
        """Make data array"""
        P, N = self.data.shape
        data_dims = [
            tiledb.Dim(name=n, domain=(1, d), tile=np.round(d/10, 0), dtype=np.int32) for \
            n, d in zip(["feature", "sample"], self.data.shape)
        ]
        data_attrs = [
            tiledb.Attr(name=self.feature, dtype=float)
        ]
        tiledb.Array.create(
            self.uris["data"].as_posix(),
            tiledb.ArraySchema(domain=tiledb.Domain(data_dims),
                               sparse=False, attrs=data_attrs)
        )

    def _make_label_arrays(self) -> None:
        label_attrs = [tiledb.Attr(name="idx", dtype=np.int32)]
        for __ in ["sample", "feature"]:
            label_dims = tiledb.Dim(name=__, tile=None, dtype="ascii")
            tiledb.Array.create(
                self.uris[__].as_posix(),
                tiledb.ArraySchema(domain=tiledb.Domain(label_dims),
                                   sparse=True, attrs=label_attrs)
            )
    
    @contextmanager
    def _open(self, mode: str = "r") -> tiledb.Array:
        with ExitStack() as stack:
            yield {k: stack.enter_context(tiledb.open(uname.as_posix(), mode)) for k, uname in self.uris.items()}
    
    def _ingest(self) -> None:
        with self._open(mode="w") as arrays:
            arrays["data"].meta["ck_version"] = self._meta()
            arrays["data"][:] = self.data.values
            arrays["sample"][self.data.columns.values] = np.arange(1, self.data.shape[1]+1)
            arrays["feature"][self.data.index.values] = np.arange(1, self.data.shape[0]+1)

    def _get(
        self,
        query_f: List[str],
        query_s: List[str],
    ) -> None:
        with self._open(mode="r") as arrays:

            fidx = arrays["feature"].df[query_f].set_index("idx")
            sidx = arrays["sample"].df[query_s].set_index("idx")
            q = arrays["data"].query(attrs=(self.feature,), coords=True)

            selected = q.multi_index[fidx.index.values, sidx.index.values]

            return (pd.DataFrame(
                selected[self.feature],
                index=selected["feature"][:,0],
                columns=selected["sample"][0,:]
            )
            .join(fidx)
            .set_index("feature")
            .transpose()
            .join(sidx)
            .set_index("sample"))

    def __getitem__(self, query):
        try:
            return self._get(query)
        except:
            raise ValueError(f"Query failed for {query}")

    @classmethod
    @contextmanager
    def read(cls, data, feature: str, tiledb_grp: Dict[str, str], cache_dir: Path = CACHE_DIR):
        return cls(data, feature, tiledb_grp, cache_dir, overwrite=False)._open(mode="r")
    

def test_get(data, db: DataDB):
    for i in range(10):
        print(i)
        sample_df = data.sample(n=10).T.sample(10).T
        qsamp = sample_df.columns[:10].tolist()
        qfeat = sample_df.index[:10].tolist()
        assert (db
                ._get(qfeat, qsamp)
                .sort_index(0).sort_index(1)
                .equals(data.loc[qfeat, qsamp].T
                        .sort_index(0).sort_index(1)))
