import os
from pathlib import Path
from typing import List, Optional
import warnings

import numpy as np
import pandas as pd
import pyensembl
import pyranges as pr

import casskit.io.utils as io_utils
import casskit.config as config
from casskit.descriptors import OneOf


class EnsemblData:
    
    # From ensembl_object.species.reference_assemblies
    # PYENSEMBL_ASSEMBLIES = {"GRCh38": 77, "GRCh37": 75, "NCBI36": 54}
    PYENSEMBL_ASSEMBLIES = {"GRCh38": 77, "GRCh37": 75}

    def __init__(
        self,
        assembly: str = "GRCh37",
        cache_dir: Optional[Path] = None
    ):
        self.assembly = assembly
        self.release = self.PYENSEMBL_ASSEMBLIES[assembly]
        if cache_dir is None:
            cache_dir = config.CACHE_DIR
        self.cache_dir = Path(cache_dir)
        
        # Set up local cache of subsetted GTF files
        os.environ['PYENSEMBL_CACHE_DIR'] = self.cache_dir.as_posix()
        if not self.gtf_path.exists():
            warnings.warn(f"Downloading GTF file for {self.assembly}...")
            self.ensembl.download()
            self.ensembl.index()

    @property
    def ensembl(self) -> pyensembl.EnsemblRelease:
        return pyensembl.EnsemblRelease(self.release)
        
    @property
    def gtf_path(self) -> Path:
        """pyensembl.EnsemblRelease.gtf_path"""
        return Path(self.cache_dir, "pyensembl", self.assembly, f"ensembl{self.release}",
                    f"Homo_sapiens.{self.assembly}.{self.release}.gtf.gz")

    @property
    def cached_subset(self) -> pr.PyRanges:
        self.set_cache(self.cache_dir / f"ensembl_{self.assembly}_{self.release}_subset.gtf")
        return self.subset_to_genes(self.gtf_path)

    @property
    def cached_tss(self) -> pr.PyRanges:
        self.set_cache(self.cache_dir / f"ensembl_{self.assembly}_{self.release}_tss.gtf")
        return self.tss_from_ensembl(pr.read_gtf(self.gtf_path))

    def set_cache(self, path_cache: Path) -> None:
        self.path_cache = path_cache
        self.read_cache = lambda cache: pr.read_gtf(cache)
        self.write_cache = lambda data, cache: data.to_gtf(cache)
        
    @io_utils.cache_on_disk
    def subset_to_genes(self, gtf_path) -> pr.PyRanges:
        """Get gene feature from Ensembl."""
        gr = pr.read_gtf(gtf_path)
        return gr[gr.Feature == "gene"]
    
    @classmethod
    def build_caches(cls):
        for assembly in cls.PYENSEMBL_ASSEMBLIES.keys():
            __ = cls(assembly).cached_subset

    @classmethod
    def to_pyranges(cls, assembly: str = "GRCh37"):
        return cls(assembly).cached_subset

    @classmethod
    def to_df(cls, assembly: str = "GRCh37"):
        return cls(assembly).cached_subset.df

    @classmethod
    def get_tss(cls, assembly: str = "GRCh37"):
        return cls(assembly).cached_tss.df

    @io_utils.cache_on_disk
    def tss_from_ensembl(self, ensembl_pr) -> pr.PyRanges:
        """Get TSSs from Ensembl.
        
        Returns the first exon of each protein-coding gene
        and pseudogene.
        """
        return (ensembl_pr[((ensembl_pr.Feature == "exon") &
                            (ensembl_pr.exon_number == "1") &
                            ((ensembl_pr.gene_biotype == "protein_coding") | 
                             (ensembl_pr.gene_biotype == "pseudogene")))]
                .df
                .astype({
                    "Chromosome": str, "Start": "int32", "End": "int32",
                    "gene_name": str, "gene_id": str
                })
                .groupby(["Chromosome", "gene_name", "gene_id"])
                .agg(
                    Start=pd.NamedAgg(column="Start", aggfunc=np.nanmin),
                    End=pd.NamedAgg(column="End", aggfunc=np.nanmax),
                )
                .reset_index()
                .pipe(pr.PyRanges, int64=True))

    @classmethod
    def annotate_pyranges(
        cls,
        gene_ids: List[str],
        assembly: str = "GRCh37",
        identifier: str = "gene_id"
    ) -> pr.PyRanges:
        OneOf("gene_id", "gene_name").validate(identifier)
        gr = cls.to_pyranges(assembly)
        return gr[getattr(gr, identifier).isin(gene_ids)]

    @classmethod
    def annotate_df(
        cls,
        gene_ids: List[str],
        assembly: str = "GRCh37",
        identifier: str = "gene_id"
    ) -> pd.DataFrame:
        return cls.annotate_pyranges(
            gene_ids, assembly=assembly, identifier=identifier
        ).df


get_ensembl_tss = EnsemblData.get_tss
""""""

build_ensembl_cache = EnsemblData.build_caches
""""""

get_ensembl = EnsemblData.to_df
""""""

annotate_genes = EnsemblData.annotate_df
""""""