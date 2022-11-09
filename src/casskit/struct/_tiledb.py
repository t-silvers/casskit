from pathlib import Path
import shutil

import numpy as np
import tiledb

cfg = tiledb.Ctx().config()
cfg.update(
  {
    'py.init_buffer_bytes': 1024**2 * 50
  }
)
tiledb.default_ctx(cfg)


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

def make_sample_dimension_labels_db(
    cache_dir: Path = Path("~/.cache").expanduser(),
    overwrite: bool = False,
):
    # Will jsut be DenseArray 'dimension labels' in the future
    uri = cache_dir / "tcga" / "samples.tiledb"

    if (Path(uri).exists() and overwrite is True):
        shutil.rmtree(uri)

    # create the schema.
    dims = [
        tiledb.Dim(name="cancer", tile=None, dtype="ascii"),
        tiledb.Dim(name="sample", tile=None, dtype="ascii"),
    ]
    attrs = [
        tiledb.Attr(name="end_pos", dtype=float, var=False),
        tiledb.Attr(name="ensembl_gene_id", dtype="S0", var=True),
        tiledb.Attr(name="value", dtype=float, var=False),
    ]
    schema = tiledb.ArraySchema(domain=tiledb.Domain(dims), attrs=attrs, sparse=True)

    # create the array.
    tiledb.SparseArray.create(uri.as_posix(), schema)

def make_feature_dimension_labels_db(
    cache_dir: Path = Path("~/.cache").expanduser(),
    overwrite: bool = False,
):
    # Will jsut be DenseArray 'dimension labels' in the future
    uri = cache_dir / "tcga" / "features.tiledb"

    if (Path(uri).exists() and overwrite is True):
        shutil.rmtree(uri)

    # create the schema.
    dims = [
        tiledb.Dim(name="chrom", tile=None, dtype="ascii"),
        tiledb.Dim(name="start_pos", domain=(1, 25E7), tile=10000, dtype=np.uint32),
        tiledb.Dim(name="feature", tile=None, dtype="ascii"),
    ]
    attrs = [
        tiledb.Dim(name="end_pos", domain=(1, 25E7), tile=10000, dtype=np.uint32),
        tiledb.Attr(name="ensembl_gene_id", dtype="S0", var=True),
    ]
    schema = tiledb.ArraySchema(domain=tiledb.Domain(dims), attrs=attrs, sparse=True)

    # create the array.
    tiledb.SparseArray.create(uri.as_posix(), schema)
