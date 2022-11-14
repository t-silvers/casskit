import latexify
import numpy as np
from numpy import linalg as LA
from sklearn.pipeline import Pipeline

from casskit.pipelines import base


class ZhangBojorquezGomezPipeline(base.CKPipeline):

    AUTHORS = [["W", "Zhang"], ["A" "Bojorquez-Gomez"], ["JF", "Kreisberg"], ["T", "Ideker"]]
    DOI = "10.1038/s41588-018-0091-2"
    URL = "https://www.nature.com/articles/s41588-018-0091-2"
    CODE = "https://github.com/wzhang1984/Noncoding-tumor-mutation-paper"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    @latexify.function
    def expression():
        # TODO: Add expression
        # cp.norm2(X @ beta - Y)**2 + lambd * cp.norm1(beta)
        pass

    @property
    def somatic_eqtl_analysis(self) -> Pipeline:
        rnaseq_pp = """
        1. RSEM-UQ
        2. RSEM-UQ >1 in >50% of tumors
        3. log2 transformed and z-score standardized
        """

        cna_pp = """
        1. GISTIC2
        2. Missing data <- 0 (diploid)
        """
        
        methylation_pp = """
        1. methylation probes mapped to the promoter regions of genes (+/-1kb from the TSS)
        2. gene <- mean(beta) of mapped probes
        """

        mutation_pp = """
        1. recurrently mutated loci
        2. promoters or enhancers <- locus +/- 100bp
        3. IF len(enhancer) > len(locus)
            gene <- intersection([promoters, enhancers], locus) > 50% of locus sequence
            OR
            gene enhancer region defined by GeneHancer
           ELSE
            gene <- intersection([promoters, enhancers], locus) > 50% of enhancers sequence
        """
        
        tumor_pp = """
        t <- tumor type, as binary variable (card(t) = 21)
        """

        ancestry_pp = """
        r <- r in {Asian; black or African American; white}, as binary variable
        """

        latent_pp = """
        1. PEER with t, r and g
        2. h determined by posterior variance of factor weights
        """

        eqtl_mapping = """
        1. solve self.expression() for beta
        2. IF somatic eQTL effect is zero (beta_1=0)
            increase lambda until beta_1 != 0
        3. pval <- F-test for nested models
        4. FDR <= 0.20, Storey's FDR correction (q-value)
        """
        
        return Pipeline(
            []
        )

    def fit(self, X, y=None):
        self.somatic_eqtl_analysis.fit(X)
        return self

    def transform(self, X):
        self.fit(X)
        return self.somatic_eqtl_analysis.transform(X)
    
    def fit_transform(self, X, y=None):
        """Custom fit_transform method for checks."""
        self.transformed = self.somatic_eqtl_analysis.fit_transform(X)
        return self.transformed