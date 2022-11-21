from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class SimTCGA:
    
    template: pd.DataFrame
    N: int = 500
    I: int = 10
    seed: int = 212
    rng: np.random.Generator = field(init=False)
    grn: pd.DataFrame = field(init=False)
    copynumber: pd.DataFrame = field(init=False)
    expression: pd.DataFrame = field(init=False)
    
    def __post_init__(self):
        super().__setattr__("rng", np.random.default_rng(self.seed))
        super().__setattr__("grn", self.make_grn())
        super().__setattr__("copynumber", self.make_copynumber())
        super().__setattr__("expression", self.make_expression())

    def make_grn(self):
        return (self.template
                .groupby("Chromosome")
                .agg({"Start": "min", "End": "max"})
                .apply(lambda x: self.rng.binomial(1, p=.4, size=self.I) * \
                       self.rng.integers(x[0], x[1], size=self.I), axis=1)
                .rename("Start")
                .apply(pd.Series)
                .melt(ignore_index=False, value_name="Start")
                .assign(gene_id=lambda x: x.pop("variable").map('egene-{}'.format),
                        beta=lambda x: self.rng.normal(0, 1, size=len(x)),)
                .reset_index()
                .query("Start != 0"))

    def make_copynumber(self):
        chroms = self.template.Chromosome.unique()
        cn_skeleton = self.chrom_swap_augmented(self.N,
                                                self.template["sample"].unique(),
                                                chroms,
                                                self.rng)

        return (cn_skeleton
                .assign(sample_sim = lambda x: [f"TCGA-00-{i:04}" for i in x.pop("sample_sim")])
                .merge(self.template)
                .drop("sample", axis=1)
                .rename(columns=dict(sample_sim="sample")))

    @staticmethod
    def chrom_swap_augmented(
        n_samples: int = 500,
        samples: list = None,
        groups: list = None,
        rng = np.random.default_rng(),
    ):
        return (
            pd.DataFrame(rng.choice(samples,
                                    replace=True,
                                    size=[len(groups), n_samples]),
                        index=groups)
            .rename_axis("Chromosome")
            .melt(value_name="sample",
                var_name="sample_sim",
                ignore_index=False)
            .reset_index()
        )

    def make_expression(self):
        design = (self.grn
                  .merge(self.copynumber, on="Chromosome", suffixes=("", "_seg"))
                  .query("Start >= Start_seg and Start <= End")
                  .pivot_table(index=["gene_id", "beta"], columns="sample", values="value", fill_value=0))

        expression = (design
                      .reset_index()
                      .groupby("gene_id")
                      .apply(lambda df: np.dot(df.beta.values, df.drop("beta", axis=1).values))
                      .apply(pd.Series))

        expression.columns = design.columns

        return expression


def simulate_tcga(
    template: pd.DataFrame,
    N = 500,
    I = 10,
    seed: int = 212,
):
    return SimTCGA(template=template, N=N, I=I, seed=seed)