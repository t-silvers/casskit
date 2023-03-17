from abc import ABC, abstractclassmethod
from typing import Dict, List

import dask.dataframe as dd
import numpy as np
import pandas as pd

from ..utils import simulate_ids


class GRNBuilder(ABC):
    def __init__(self, data):
        self.data = data

    @abstractclassmethod
    def simulate_grn_eqtls(self, **kwargs):
        pass

    @abstractclassmethod
    def simulate_eqtl_design(self, **kwargs):
        pass

class cneQTLGRNBuilder(GRNBuilder):
    def __init__(self, data):
        super().__init__(data)

    def simulate_grn_eqtls(
        self,
        n_genes,
        *args,
        avg_num_cneqtls=12,
        max_per_chrom=3,
        random_state=None,
        **kwargs,
    ):
        rng = np.random.default_rng(random_state)
        p_cneqtl = avg_num_cneqtls / max_per_chrom / self.data["Chromosome"].nunique()
        
        # TODO: Change input requirements to factor out initial parsing.
        #       Pass data and possible eqtls separately?
        cneqtls = (self.data
            .groupby("Chromosome")
            .agg({"Start": "min", "End": "max"})
            .apply(lambda x: rng.binomial(1, p=p_cneqtl, size=(max_per_chrom, n_genes)) * \
                rng.integers(x[0], x[1], size=(max_per_chrom, n_genes)), axis=1)
            .explode()
            .apply(pd.Series)
            .melt(ignore_index=False, value_name="Start")
            .assign(gene_id=lambda x: simulate_ids(x.pop("variable"), "gene_"))
            .query("Start != 0")
            .set_index("gene_id", append=True)
            .add_suffix("_cneqtl")
            .reset_index()
            .assign(eqtl=lambda x: \
                x["Chromosome"] + "_" + x["Start_cneqtl"].astype(str)))

        self.eqtls = cneqtls

    def simulate_eqtl_design(self, *args, **kwargs):
        # Required for dask pivoting
        data = self.data.copy()
        eqtls = self.eqtls.copy()
        data["sample_id"] = data["sample_id"].astype("category")
        eqtls["gene_id__eqtl__Chromosome__Position"] = \
            (eqtls["gene_id"] + "__" 
             + eqtls["eqtl"] + "__" 
             + eqtls["Chromosome"] + "__" 
             + eqtls["Start_cneqtl"].astype(str)
             ).astype("category")

        # Convert to dask
        cneqtls_dd = dd.from_pandas(eqtls, npartitions=4)
        copynumber_dd = dd.from_pandas(data, npartitions=4)
        cneqtl_design_dd = dd.merge(cneqtls_dd, copynumber_dd, on="Chromosome")
        cneqtl_design_dd = (cneqtl_design_dd
                        .query("Start_cneqtl >= Start & Start_cneqtl <= End")
                        .pivot_table(index="gene_id__eqtl__Chromosome__Position",
                                    columns="sample_id",
                                    values="value"))
        cneqtl_design = cneqtl_design_dd.compute()

        # Fix index from dask constraints
        predask_idx = cneqtl_design.index.str.split("__", expand=True)
        predask_idx.names = ["gene_id", "eqtl", "Chromosome", "Position"]
        
        # Restore dtypes for numeric
        position_idx = predask_idx.levels[3].astype(int)
        predask_idx = predask_idx.set_levels(position_idx, level="Position")
        
        # Set index
        cneqtl_design = cneqtl_design.set_index(predask_idx)
        self.design = cneqtl_design

class muteQTLGRNBuilder(GRNBuilder):
    def __init__(self, data):
        super().__init__(data)

    def simulate_grn_eqtls(self, n_genes, *args, avg_num_muteqtls=2, random_state=None, **kwargs):
        rng = np.random.default_rng(random_state)
        possible_muteqtls = self.data.columns
        p_cneqtl = avg_num_muteqtls / possible_muteqtls.size
        if p_cneqtl > 0:
            selected_muteqtls = rng.binomial(1, p=p_cneqtl, size=(possible_muteqtls.size, n_genes))

            muteqtls = []
            for i, gene in enumerate(simulate_ids(n_genes, "gene_")):
                gene_muteqtls = possible_muteqtls[selected_muteqtls[:, i].astype(bool)]
                gene_muteqtls = [(gene, v) for v in gene_muteqtls]
                muteqtls = muteqtls + gene_muteqtls
            
            muteqtls = pd.DataFrame.from_records(muteqtls, columns=["gene_id", "eqtl"])
            self.eqtls = muteqtls
        
        else:
            self.eqtls = pd.DataFrame(columns=["gene_id", "eqtl"])

    def simulate_eqtl_design(self, *args, **kwargs):
        muteqtl_design = (self.data
                    .melt(var_name="eqtl", ignore_index=False)
                    .query("value != 0")
                    .rename_axis("sample_id")
                    .reset_index()
                    .merge(self.eqtls)
                    # TEMP
                    .assign(Position=0, Chromosome="chr1")
                    .pivot_table(index=["gene_id", "eqtl", "Chromosome", "Position"],
                                columns="sample_id",
                                values="value",
                                fill_value=0))

        self.design = muteqtl_design

class GRNgineer:    
    def __init__(self, min_var=0.1, min_herit=None, unif_herit=None, **kwargs):
        self._design = None
        self._grn = None
        self.min_var = min_var
        self.min_herit = min_herit
        self.unif_herit = unif_herit

    @property
    def design(self):
        return self._design

    @design.setter
    def design(self, design):
        # Filter for low variance
        locus_var = design.var(axis=1)
        passing_loci = locus_var.mask(lambda s: s < self.min_var)
        design_filt = design.loc[passing_loci.dropna().index, :]
        design_all = pd.concat([self._design, design_filt], axis=0)
        # If some samples not in new design, fill with 0
        design_all = design_all.fillna(0)
        self._design = design_all.sort_index()

    @property
    def grn(self):
        return self._grn

    @grn.setter
    def grn(self, random_state):
        rng = np.random.default_rng(random_state)
        grn = self.design.index.to_frame(index=True).filter([])
        if self.min_herit is not None:
            betas = simulate_eqtl_effects(self.design.shape[0], random_state)
            herit = (self.design.var(axis=1) * betas**2).values
            grn = grn.assign(beta=betas, herit=herit)
            # Filter by minimum herit
            grn = grn.query("herit > @self.min_herit")
        if self.unif_herit is not None:
            herit = self.unif_herit
            betas = np.sqrt(herit / self.design.var(axis=1))
            betas_sign = rng.choice([-1, 1], size=betas.size)
            betas = np.multiply(betas, betas_sign)
        else:
            betas = simulate_eqtl_effects(self.design.shape[0], random_state)
            herit = (self.design.var(axis=1) * betas**2).values
            pass
        grn = grn.assign(beta=betas, herit=herit)
        self._grn = grn

    def simulate_grn(
        self,
        builders: Dict[str, GRNBuilder],
        n_genes,
        *args,
        random_state=None,
        **kwargs,
    ):
        for mol, builder in builders.items():
            builder.simulate_grn_eqtls(n_genes, **kwargs)
            builder.simulate_eqtl_design(**kwargs)
            design = (builder.design
                      .assign(mol=mol)
                      .set_index("mol", append=True)
                      .swaplevel("eqtl", "mol"))
            self.design = design
        self.grn = random_state

    @classmethod
    def from_datasets(
        cls,
        datasets: Dict[str, pd.DataFrame],
        n_genes: int,
        *args,
        **kwargs
    ):
        builders = dict(copynumber=cneQTLGRNBuilder(datasets.get("copynumber", None)),
                        mutation=muteQTLGRNBuilder(datasets.get("mutation", None)))
        builders = {k: v for k, v in builders.items() if v.data is not None}
        grn_engineer = cls(**kwargs)
        grn_engineer.simulate_grn(builders=builders, n_genes=n_genes, **kwargs)
                
        return grn_engineer.grn, grn_engineer.design
    
def simulate_eqtl_effects(size, random_state):
    rng = np.random.default_rng(random_state)
    # TODO: Add cis-cn-eQTL
    # TODO: Other dists
    # TODO: Refactor to include in builder classes
    betas = rng.normal(0, 1, size)
    return betas
        