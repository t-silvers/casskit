from dataclasses import dataclass, field
import warnings

import numpy as np
import pandas as pd

# TODO: Refactor with builder class

@dataclass(frozen=True)
class SimTCGA:
    template: pd.DataFrame
    N: int = 500
    I: int = 10
    seed: int = 212
    augmented: bool = False
    min_h2: float = .1
    rng: np.random.Generator = field(init=False)
    grn: pd.DataFrame = field(init=False)
    copynumber: pd.DataFrame = field(init=False)
    expression: pd.DataFrame = field(init=False)
    
    def __post_init__(self):
        super().__setattr__("rng", np.random.default_rng(self.seed))

        # Make copy number
        # ----------------
        super().__setattr__("copynumber", self.make_copynumber())
        
        # Make GRN
        # --------
        proposal_grn = self.make_grn()

        # Filter for identifiability
        accepted_grn = self._filter_identifiability(proposal_grn)
        super().__setattr__("grn", accepted_grn)
        
        # Make expression
        # ---------------
        super().__setattr__("expression", self.make_expression())

    def make_copynumber(self):
        if self.augmented is True:
            chroms = self.template.Chromosome.unique()
            cn_skeleton = chrom_swap_augmented(
                self.N,
                self.template["sample"].unique(),
                chroms,
                self.rng
            )

            return (cn_skeleton
                    .assign(sample_sim = lambda x: [f"TCGA-00-{i:04}" for i in x.pop("sample_sim")])
                    .merge(self.template)
                    .drop("sample", axis=1)
                    .rename(columns=dict(sample_sim="sample")))
        
        else:
            n_samples = self.template["sample"].unique()
            n_sim = self.N
            if len(n_samples) < self.N:
                warnings.warn(f"Requested {self.N} samples, "
                              f"but only {len(n_samples)} are available.")
                n_sim = len(n_samples)

            samples = self.rng.choice(self.template["sample"].unique(),
                                      size=n_sim,
                                      replace=False)
            
            return self.template.query("sample in @samples")

    def make_grn(self, per_chrom=True, cis=True):
        if per_chrom is True:
            # Sampling twice is a hacky way to have multiple cn-eQTLs on the
            # same chromosome with Pr=p**2=0.25
            cneqtls1 = (self.template
                       .groupby("Chromosome")
                       .agg({"Start": "min", "End": "max"})
                       .apply(lambda x: self.rng.binomial(1, p=.5, size=self.I) * \
                           self.rng.integers(x[0], x[1], size=self.I), axis=1))

            cneqtls2 = (self.template
                       .groupby("Chromosome")
                       .agg({"Start": "min", "End": "max"})
                       .apply(lambda x: self.rng.binomial(1, p=.5, size=self.I) * \
                           self.rng.integers(x[0], x[1], size=self.I), axis=1))

            # combine
            cneqtls = pd.concat([cneqtls1, cneqtls2], axis=0)
        
        else:
            raise NotImplementedError("Not implemented yet.")

        if cis is True:
            # make cis, cis = abs(10 * trans)
            cis_cneqtl = cneqtls.sample(1).abs().mul(10)
            cneqtls = pd.concat([cis_cneqtl, cneqtls.drop(cis_cneqtl.index)], axis=0)

        return (cneqtls
                .rename("Start")
                .apply(pd.Series)
                .melt(ignore_index=False, value_name="Start")
                    
                # Simulate cn-eQTL effects from normal distribution
                .assign(gene_id=lambda x: x.pop("variable").map('egene-{}'.format),
                        beta=lambda x: self.rng.normal(0, 0.5, size=len(x)),)
                .reset_index()
                .query("Start != 0"))

    def _filter_identifiability(self, grn) -> pd.DataFrame:
        """Filter GRN for identifiability based on ~h2."""
        proposal_design = make_design(grn, self.copynumber)
        proposal_betasq = proposal_design.index.get_level_values("beta")**2
        filtered_grn = (proposal_design
                        .assign(hfilt=\
                            lambda x: x.var(axis=1) \
                                * proposal_betasq >= self.min_h2)
                        .query("hfilt")
                        .drop("hfilt", axis=1)
                        .index.to_frame(index=False)
                        .merge(grn))
        
        # Record h2
        filtered_design = make_design(filtered_grn, self.copynumber)
        filtered_betasq = filtered_design.index.get_level_values("beta")**2
        filtered_grn["herit"] = (filtered_design.var(axis=1)
                                 * filtered_betasq).values
        
        return filtered_grn

    def make_expression(self):
        design = make_design(self.grn, self.copynumber)
        expression = (design
                      .droplevel(["Chromosome", "herit"])
                      .reset_index()
                      .groupby("gene_id")
                      .apply(lambda df: np.dot(df.beta.values, df.drop("beta", axis=1).values))
                      .apply(pd.Series))

        # Add noise
        rng = np.random.default_rng()
        esig = rng.normal(0, 1, expression.shape)
        expression += esig
        
        expression.columns = design.columns

        return expression.T

    def __repr__(self):
        return (f"SimTCGA(N={self.N}, I={self.I}, "
                "augmented={self.augmented}, seed={self.seed})")

def make_design(grn, copynumber):
    # idx_cols = ["gene_id", "Chromosome", "herit", "beta"]
    idx_cols = [col for col in grn.columns if 
                col not in ["Start", "End", "sample", "value"]]
    return (grn
            .merge(copynumber,
                    on="Chromosome",
                    suffixes=("", "_seg"))
            .query("Start >= Start_seg and Start <= End")
            .pivot_table(index=idx_cols,
                         columns="sample",
                         values="value",
                         fill_value=0))

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

def simulate_tcga(
    template: pd.DataFrame,
    N = 500,
    I = 10,
    seed: int = 212,
    augmented: bool = False,
    min_h2: int = .1,
):
    return SimTCGA(template=template,
                   N=N,
                   I=I,
                   seed=seed,
                   augmented=augmented,
                   min_h2=min_h2)