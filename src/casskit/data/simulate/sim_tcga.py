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


class SimTCGA:
    def __init__(
        self,
        cancer,
        I: int = 100,
        N: int = 100,
        P: int = 500,
        seed: int = None
    ):
        self.cancer = cancer
        self.I = I
        self.N = N
        self.P = P
        self.seed = seed
        self.cneqtl_size = 1E5
        self.cn_method = "swap_augment"
        self.k_cis = 2
        self.k_var = 5
        self.k_trans = 10

        self.tcga_expression = get_tcga("htseq_counts", cancer)
        self.tcga_cn = get_tcga("cnv", cancer)

    @property
    def egenes(self) -> Dict:
        return self.get_candidate_egenes()

    @property
    def chroms(self) -> List:
        return self.tcga_cn.Chrom.unique().tolist()

    @property
    def candidate_variants(self) -> pd.DataFrame:
        return self.variants.columns.tolist()

    @property
    def candidate_cneqtls(self) -> Dict:
        return self.get_candidate_cneqtls()

    @property
    def variants(self) -> pd.DataFrame:
        return simulate_variants(self.N, self.P)

    @property
    def copynumber(self) -> pd.DataFrame:
        return simulate_copynumber(
            self.N, self.P, cn_method=self.cn_method, copynumber_template=self.tcga_cn
        )

    @property
    def expression(self) -> pd.DataFrame:
        grn, expression = self.make_expression()
        self.grn = grn
        return expression

    def get_candidate_egenes(self) -> np.ndarray:
        egenes = self.tcga_expression.Ensembl_ID.str.split(".", expand=True)[0].unique()
        egene_tss = get_ensembl_tss()
        
        return (egene_tss
                .query("gene_id in @egenes")
                .sample(n=self.I, random_state=self.seed, replace=False)
                .filter(["gene_id", "gene_name", "Chromosome", "Start", "End"])
                .to_dict("records"))

    def get_candidate_cneqtls(self, coarsen=-4) -> Dict:
        tcga_cn_start = self.tcga_cn.Start.round(coarsen).astype(int)
        cn_bin_starts = np.arange(tcga_cn_start.min(), tcga_cn_start.max()+self.cneqtl_size, self.cneqtl_size)

        return (self.tcga_cn
                .assign(
                    cn_bin_start=pd.cut(
                        tcga_cn_start, bins=cn_bin_starts, labels=cn_bin_starts[:-1]
                    ).astype(float),
                ).drop_duplicates(subset=["Chrom", "cn_bin_start"])
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
            regulators=grn_sim,
            variants=self.variants,
            copynumber=self.copynumber
        )

    def make_expression(self):
        grns, expressions = [], []
        for egene in tqdm.tqdm(self.egenes):
            print(egene)
            grn, grn_df = self.simulate_egene_grn(
                egene.get("gene_id"), egene.get("gene_name"),
                (egene.get("Chromosome"), egene.get("Start"), egene.get("End"))
            )
            expression = self.simulate_egene_expression(
                egene.get("gene_id"), grn
            )
            grns.append(grn_df)
            expressions.append(expression)

        return pd.concat(grns, axis=1), pd.concat(expressions, axis=1)

simulate_tcga = lambda cancer: SimTCGA(cancer, I=100, N=100, P=500, seed=1234)