from casskit.data.simulate.sim_variants import simulate_variants

class SimTCGA:
    
    def __init__(
        self,
        egenes,
        cancer,
        I=100,
        N=100,
        P=500,
        template_cn: pd.DataFrame = None,
        cn_method="swap_augment",
        k_cis: int = 0,
        k_var: int = 2,
        k_trans: int = 10,
        cneqtl_size: int = int(1E5),
        seed: int = 212,
    ):
        self.cancer = cancer
        self.egenes = egenes
        self.I = I
        self.N = N
        self.P = P
        self.seed = seed
        self.template_cn = template_cn
        self.k_cis = k_cis
        self.k_var = k_var
        self.k_trans = k_trans
        self.cneqtl_size = cneqtl_size

        self.chroms = template_cn.index.get_level_values("Chromosome").unique().tolist()
        self.candidate_cneqtls = self.cneqtl_candidates(template_cn)

        if cn_method == "original":
            self.copynumber = template_cn

        elif cn_method == "swap_augment":
            samples = template_cn.columns.tolist()

            return self.chrom_swap_augmented(
                template_cn, self.N, samples, self.chroms
            )

        self.variants = simulate_variants(self.N, self.P)
        self.candidate_variants = self.variants.columns.tolist()


    def get_candidate_egenes(self) -> np.ndarray:
        egene_tss = get_ensembl_tss().df
        
        return (egene_tss
                .query("gene_id in @self.egenes & Chromosome in @self.chroms")
                .sample(n=self.I, random_state=self.seed, replace=False)
                .filter(["gene_id", "gene_name", "Chromosome", "Start", "End"])
                .to_dict("records"))
            
    def chrom_swap_augmented(
        self,
        copynumber: pd.DataFrame,
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
            .rename_axis("Chrom")
            .melt(value_name="sample",
                var_name="sample_sim",
                ignore_index=False)
            .reset_index()
            .merge(copynumber)
            .drop("sample", axis=1)
            .rename(columns={"sample_sim": "sample"})
        )

    def cneqtl_candidates(self, data):
        """Get candidate loci for cn-eQTL simulations.
        
        Filter by high-variance and space between loci.
        """
        return (data.var(axis=1)
                .rename("var")
                .reset_index()
                .reset_index()
                .mask(lambda x: x["var"] < x["var"].mean())
                .dropna()
                .groupby("Chromosome", group_keys=False)
                .apply(lambda df: df.assign(ix_diff = lambda x: x["index"].rolling(2).apply(np.diff)).fillna(100))
                .query("ix_diff > 5")
                .filter(["Chromosome", "Start", "End"]))


