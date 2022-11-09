from ast import Expression
from contextlib import ExitStack
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
    'py.init_buffer_bytes': 1024**2 * 500
  }
)
tiledb.default_ctx(cfg)


def _mkdim(name: str, n: int = 1E5, t: int = 100) -> tiledb.Dim:
    return tiledb.Dim(name=name, domain=(1, n), tile=t, dtype=np.int32)

def _mkattr(name: str) -> tiledb.Attr:
    return tiledb.Attr(name=name, dtype=float)

def _uri(cache_dir: Path, groups: Dict[str, str], dbtype: str,
         prefix: Union[str, Path] = "", **params) -> Path:
    
    prefix = Path(prefix)
    for k,v in groups.items():
        prefix /= f"{k}~{v}"
        try:
            tiledb.group_create(f"{k}~{v}")
        except TileDBError as err:
            print("TileDB exception: ", err)
        
    prefix /= ".".join(f"{k}~{v}" for k,v in params.items())
    uri = cache_dir / f"{prefix}.{dbtype}.tiledb"
    uri.parent.mkdir(parents=True, exist_ok=True)
    
    return uri



uri = _uri(ck.CACHE_DIR, {"cancer": cancer}, mol="expression")
dims = [_mkdim("feature"), _mkdim("sample")]
attrs = [_mkattr(feature)]

data2db(uri, dims, attrs)



# should first make label description db
# then make data db


def data2db(
    data,
    uri: Path,
    cache_dir: Path = CACHE_DIR,
    overwrite: bool = False,
    **kwargs
) -> None:
    
    
    cancer = "TCGA-CHOL"
    feature = "htseq_counts"
    # cancer = "GDC-PANCAN"
    
    data = ck.io.get_tcga(feature, cancer).set_index("Ensembl_ID")
    
    cache_dir: Path = ck.CACHE_DIR
    data_attrs = [tiledb.Attr(name=feature, dtype=float)]
    kwargs = {"mol": "expression"}
    groups = {"cancer": cancer}
    n = 1E5
    t = 10000
    chunksize = 10000
    source = "UCSC Xena"
    dom_start = 1
    overwrite = True
    
    
    # Make URIs
    uris = {}
    for __ in ["data", "sample", "feature"]:
        uris[__] = _uri(cache_dir, groups, dbtype=__, **kwargs)
        if (uris[__].exists() and overwrite is True):
            shutil.rmtree(uris[__])

    # Make data array
    data_dims = [
        tiledb.Dim(name="sample", domain=(dom_start, n), tile=t, dtype=np.int32),
        tiledb.Dim(name="feature", domain=(dom_start, n), tile=t, dtype=np.int32)
    ]

    tiledb.Array.create(uris["data"].as_posix(),
                        tiledb.ArraySchema(domain=tiledb.Domain(data_dims),
                                           sparse=False, attrs=data_attrs))

    # Make label description arrays
    label_attrs = [tiledb.Attr(name="idx", dtype=np.int32)]
    for __ in ["sample", "feature"]:
        label_dims = tiledb.Dim(name=__, tile=None, dtype="ascii")
        tiledb.Array.create(uris[__].as_posix(),
                            tiledb.ArraySchema(domain=tiledb.Domain(label_dims),
                                               sparse=True, attrs=label_attrs))

    # Ingest data
    with ExitStack() as stack:
        dbs = {k: stack.enter_context(tiledb.open(uname.as_posix(), 'w')) for k, uname in uris.items()}

        dbs["data"].meta["source"] = source
        dbs["data"].meta["tdb_version"] = f"TileDB v{tiledb.__version__}"
        dbs["data"].meta["ck_version"] = f"CAss-Kit v{__version__}"

        for ix, chunk in data.groupby(np.arange(len(data)) // chunksize):
            print(f"Writing chunk {ix} ...")
            P, N = chunk.shape
            ix0_s, ix0_e = P*ix+dom_start, P*(ix+dom_start)+dom_start
            ix1_s, ix1_e = dom_start, (N+dom_start)
            
            dbs["data"][ix0_s:ix0_e, ix1_s:ix1_e] = chunk.values
            dbs["feature"][chunk.index.values] = np.arange(ix0_s, ix0_e)
            if ix == 0:
                dbs["sample"][chunk.columns.values] = np.arange(ix1_s, ix1_e)




with tiledb.open('/scratch/users/tsilvers/ceqtl_selection/.cache/casskit/cancer~TCGA-CHOL/mol~expression.data.tiledb', "r") as A:
    print(A.df[:])


test_samples = ["TCGA-W5-AA34-11A", "TCGA-W5-AA2U-11A"]
test_features = ["ENSG00000281131.1"]

with ExitStack() as stack:
    dbs = {k: stack.enter_context(tiledb.open(uname.as_posix(), mode="r")) for k, uname in uris.items()}
    feature_idx = dbs["feature"].df[test_features].idx.tolist()
    sample_idx = dbs["sample"].df[test_samples].idx.tolist()
    # test_data = dbs["data"][feature_idx, sample_idx]
    # test_data = dbs["data"][0:2, 1:3]
    print(dbs["data"].df[:])

with tiledb.open(uri.as_posix(), 'r') as A:
    print(A.df[:])


# The order of the dimensions in the array schema matters.
# More selective dimensions (i.e., with greater pruning power)
# should be defined before less selective ones.


with tiledb.open(uris["sample"].as_posix(), 'r') as A:
    print(A.df[[test_sample, "TCGA-W5-AA2U-11A"]])


def make_feature_labels_db(uri: Path, overwrite: bool = False,) -> None:
    
    if (Path(uri).exists() and overwrite is True):
        shutil.rmtree(uri)

    tiledb.Array.create(
        Path(uri).as_posix(),
        tiledb.ArraySchema(
            domain=tiledb.Domain([tiledb.Dim(name="feature_id", tile=None, dtype="ascii"),]),
            sparse=True,
            attrs=[tiledb.Attr(name="feature", dtype="ascii"),]
        )
    )




@dataclass
class DataDB:
    data: pd.DataFrame

    def _get(self, **params):
        return self.db.get(**params)

    def __getitem__(self, **params):
        maf = self._get(**params)
        if maf:
            return self._type(maf)
        else:
            raise KeyError('This variant is not in the db')

    def __contains__(self, **params):
        return self._get(**params) is not None

    def __post_init__(self):
        self.db = tiledb.SparseArray(self.uri.as_posix(), mode='w')


    



def make_data_db(uri: Path, overwrite: bool = False, n_features: int = 1E5,
                 n_samples: int = 1E5,) -> None:
    
    if (Path(uri).exists() and overwrite is True):
        shutil.rmtree(uri)

    tiledb.Array.create(
        Path(uri).as_posix(),
        tiledb.ArraySchema(
            domain=tiledb.Domain(
                [
                    tiledb.Dim(name="feature", domain=(1, n_features), tile=1000, dtype=np.int32),
                    tiledb.Dim(name="sample", domain=(1, n_samples), tile=100, dtype=np.int32),
                ]
            ),
            sparse=False,
            attrs=[tiledb.Attr(name="original", dtype=float),]
        )
    )


# ==============================================================================

class VariantDB:

    def __init__(self, path, rocksdb_options=None):
        if rocksdb_options is None:
            rocksdb_options = rocksdb.Options(
                create_if_missing=True,
                max_open_files=100,
            )

        self.db = rocksdb.DB(
            path,
            rocksdb_options,
            read_only=True
        )

    @staticmethod
    def _variant_to_byte(variant):
        return bytes(str(variant), 'utf-8')

    def _type(self, value):
        raise NotImplementedError()

    def _get(self, variant):
        if not variant.startswith('chr'):
            variant = 'chr%s' % variant
        return self.db.get(self._variant_to_byte(variant))

    def __getitem__(self, variant):
        maf = self._get(variant)
        if maf:
            return self._type(maf)
        else:
            raise KeyError('This variant is not in the db')

    def __contains__(self, variant):
        return self._get(variant) is not None

    def get(self, variant, default=None):
        try:
            return self[variant]
        except KeyError:
            return default


class VariantMafDB(VariantDB):

    def _type(self, value):
        return float(value)
    