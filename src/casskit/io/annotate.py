from collections import namedtuple
import logging
import os
from pathlib import Path
from typing import List, Optional
import warnings

import pandas as pd
import pyensembl
import pyranges as pr

import casskit.io.utils as io_utils
import casskit.config as config

cache_dir = config.CACHE_DIR


class EnsemblData:
    
    # ensembl_object.species.reference_assemblies
    PYENSEMBL_ASSEMBLIES = {"GRCh38": 77, "GRCh37": 75, "NCBI36": 54}
    
    def __init__(
        self,
        assembly: str = "GRCh37",
        release: Optional[int] = None,
        cache_dir: Optional[Path] = cache_dir
    ):
        self.assembly = assembly
        if release is None:
            self.release = self.PYENSEMBL_ASSEMBLIES[assembly]
        self.set_cache(cache_dir)

    @property
    def ensembl(self) -> pyensembl.EnsemblRelease:
        ensembl = pyensembl.EnsemblRelease(self.release)
        if not ensembl.required_local_files_exist():
            warnings.warn(f"Downloading GTF file for {self.assembly}...")
            ensembl.download()
            ensembl.index()
        
        return ensembl

    @property
    def gtf_path(self) -> Path:
        """Path to GTF file.
        
        pyensembl.EnsemblRelease.gtf_path not working ...
        """
        return Path(self.cache_dir, "pyensembl", self.assembly, f"ensembl{self.release}",
                    f"Homo_sapiens.{self.assembly}.{self.release}.gtf.gz")

    @property
    def cached_subset(self) -> pr.PyRanges:
        return self.subset_to_genes()
        
    def set_cache(self, cache_dir: Path) -> None:
        self.cache_dir = Path(cache_dir)
        
        # Set pyensembl cache directory
        os.environ['PYENSEMBL_CACHE_DIR'] = self.cache_dir.as_posix()
        
        # Set up local cache of subsetted GTF files
        self.path_cache = self.cache_dir / f"ensembl_{self.assembly}_{self.release}_subset.gtf"
        self.read_cache = lambda cache: pr.read_gtf(cache)
        self.write_cache = lambda data, cache: data.to_gtf(cache)
        
    @io_utils.cache_on_disk
    def subset_to_genes(self) -> pr.PyRanges:
        """Get gene feature from Ensembl.
        """
        gr = pr.read_gtf(self.gtf_path)
        return gr[gr.Feature == "gene"]
    
    @classmethod
    def build_caches(cls):
        for assembly, release in cls.PYENSEMBL_ASSEMBLIES.items():
            __ = cls(assembly, release).cached_subset

    @classmethod
    def to_pyranges(cls, assembly: str):
        return cls(assembly).cached_subset

    @classmethod
    def to_df(cls, assembly: str):
        return cls(assembly).cached_subset.df

    @classmethod
    def annotate_pyranges(cls, assembly: str, identifier:str, gene_ids: List[str]) -> pr.PyRanges:
        gr = cls.to_pyranges(assembly)
        if identifier == "ensembl_gene_id":
            return gr[gr.gene_id.isin(gene_ids)]

        elif identifier == "hgnc_symbol":
            return gr[gr.gene_name.isin(gene_ids)]

    @classmethod
    def annotate_df(cls, assembly: str, identifier:str, gene_ids: List[str]) -> pd.DataFrame:
        return cls.annotate_pyranges(assembly, identifier, gene_ids).df


build_ensembl_cache = EnsemblData.build_caches
""""""

get_ensembl = EnsemblData.to_df
""""""

annotate_genes = EnsemblData.annotate_df
""""""

# Need resource for TSSs
