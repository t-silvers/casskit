from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import tqdm

from casskit.io import get_tcga, get_ensembl_tss
from casskit.data.simulate.sim_variants import simulate_variants
from casskit.data.simulate.sim_copynumber import simulate_copynumber
from casskit.data.simulate.sim_grn import simulate_grn
from casskit.data.simulate.sim_expression import simulate_expression
from casskit.preprocess.units import ToCounts


@dataclass(frozen=True)
class SimTCGA:

    cancer: str
    I: int
    N: int
    P: int
    seed: int = 212

    cneqtl_size: int = int(1E5)
    cn_method: str = "swap_augment"
    k_cis: int = 0
    k_var: int = 2
    k_trans: int = 10

    tcga_expression: pd.DataFrame = field(init=False, default=None)
    tcga_cn: pd.DataFrame = field(init=False, default=None)

    def __post_init__(self):
        if self.tcga_expression is None:
            super().__setattr__("tcga_expression", get_tcga("htseq_counts", self.cancer))
        
        if self.tcga_cn is None:
            super().__setattr__(
                "tcga_cn",
                (get_tcga("cnv", self.cancer)
                 .assign(
                     value=lambda x: ToCounts("log2(copy-number/2)").fit_transform(x.pop("value"))
                 ))
            )
        
        super().__setattr__("chroms", self.tcga_cn.Chrom.unique().tolist())
        super().__setattr__("egenes", self.get_candidate_egenes())
        
        # Simulate omics
        super().__setattr__("candidate_cneqtls", self.get_candidate_cneqtls())
        super().__setattr__("variants", simulate_variants(self.N, self.P))
        super().__setattr__("candidate_variants", self.variants.columns.tolist())
        super().__setattr__("copynumber",
            simulate_copynumber(
                self.N, self.P, cn_method=self.cn_method,
                copynumber_template=self.tcga_cn
            )
        )
        
        # Simulate expression
        grns, grn_df = self.make_grns()
        super().__setattr__("grns", grns)
        super().__setattr__("grn_df", grn_df)
        super().__setattr__("expression", self.make_expression(grns))
    
    def get_candidate_egenes(self) -> np.ndarray:
        egenes = self.tcga_expression.Ensembl_ID.str.split(".", expand=True)[0].unique()
        egene_tss = get_ensembl_tss().df
        
        return (egene_tss
                .query("gene_id in @egenes & Chromosome in @self.chroms")
                .sample(n=self.I, random_state=self.seed, replace=False)
                .filter(["gene_id", "gene_name", "Chromosome", "Start", "End"])
                .to_dict("records"))

    def get_candidate_cneqtls(self, coarsen=-4) -> Dict:
        tcga_cn_start = self.tcga_cn.Start.round(coarsen).astype(int)
        cn_bin_starts = np.arange(tcga_cn_start.min(), tcga_cn_start.max()+self.cneqtl_size, self.cneqtl_size)

        # For filtering by high variance cns
        N = self.tcga_cn["sample"].nunique()
        N_cutoff = int(N/10)
        var_cutoff = 1

        return (self.tcga_cn
                .assign(
                    cn_bin_start=pd.cut(
                        tcga_cn_start, bins=cn_bin_starts, labels=cn_bin_starts[:-1]
                    ).astype(float),
                )
                # filter by high variance cns
                .groupby(["Chrom", "cn_bin_start"])
                [["value"]].agg({"value": ["var", "count"]})
                ["value"]
                .query("count > @N_cutoff & var > @var_cutoff")
                .index.to_frame(index=False)
                .groupby("Chrom", group_keys=False)
                ["cn_bin_start"]
                .apply(lambda x: x.tolist())
                .to_dict())

    def simulate_egene_grn(self, gene_id, gene_name, gene_coords):
        return simulate_grn(
            gene_id,
            gene_name,
            egene_coords=gene_coords,
            chroms=self.chroms,
            var_ids=self.candidate_variants,
            copynumber_ids=self.candidate_cneqtls,
            k_cis=self.k_cis,
            k_var=self.k_var,
            k_trans=self.k_trans,
            return_df=True
        )

    def simulate_egene_expression(self, gene_id, grn_sim):
        return simulate_expression(
            gene_id,
            N=self.N,
            p=self.P,
            regulators=grn_sim,
            variants=self.variants,
            copynumber=self.copynumber
        )

    def make_grns(self):
        grns, grn_dfs = {}, []
        for egene in self.egenes:
            grn, grn_df = self.simulate_egene_grn(
                egene.get("gene_id"), egene.get("gene_name"),
                (egene.get("Chromosome"), egene.get("Start"), egene.get("End"))
            )
            grns[egene.get("gene_id")] = grn
            grn_dfs.append(grn_df.assign(gene_id=egene.get("gene_id")))

        return grns, pd.concat(grn_dfs)

    def make_expression(self, grns):
        return pd.concat(
            [self.simulate_egene_expression(
                egene.get("gene_id"), grns[egene.get("gene_id")]
             ) for egene in tqdm.tqdm(self.egenes)],
            axis=1
        )

    def __repr__(self) -> str:
        return f"SimTCGA(cancer={self.cancer}, I={self.I}, N={self.N}, P={self.P}, seed={self.seed})"

def simulate_tcga(cancer, I=100, N=100, P=500, seed=212):
    return SimTCGA(cancer, I, N, P, seed=seed)