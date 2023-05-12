from ..io import get_tcga
from ..omic import (
    CopyNumberVariation,
    MessengerRNA,
    SomaticMutation,
)
from ..multiomic import BaseMultiOmic


class TCGADataSet(BaseMultiOmic):
    """Light preprocessing."""
    def __init__(self, cancer):
        self.cancer = cancer

    ######################
    ## COPY NUMBER DATA ##
    ######################

    @property
    def cnv_data(self):
        prepared_cnv = prepare_cnv(
            get_tcga("cnv", self.cancer)
        )
        cnv_data = CopyNumberVariation(prepared_cnv)
        return cnv_data

    @property
    def gistic_data(self):
        prepared_gistic = prepare_gistic(
            get_tcga("gistic", self.cancer)
        )
        gistic_data = CopyNumberVariation(prepared_gistic)
        return gistic_data

    @property
    def masked_cnv_data(self):
        prepared_masked_cnv = prepare_masked_cnv(
            get_tcga("masked_cnv", self.cancer)
        )
        masked_cnv_data = CopyNumberVariation(prepared_masked_cnv)
        return masked_cnv_data

    #####################
    ## EXPRESSION DATA ##
    #####################

    @property
    def htseq_counts_data(self):
        prepared_htseq_counts = prepare_htseq_counts(
            get_tcga("htseq_counts", self.cancer)
        )
        htseq_counts_data = MessengerRNA(prepared_htseq_counts)
        return htseq_counts_data

    @property
    def htseq_fpkm_uq_data(self):
        prepared_htseq_fpkm_uq = prepare_htseq_fpkm_uq(
            get_tcga("htseq_fpkm-uq", self.cancer)
        )
        htseq_fpkm_uq_data = MessengerRNA(prepared_htseq_fpkm_uq)
        return htseq_fpkm_uq_data

    ###################
    ## MUTATION DATA ##
    ###################

    @property
    def muse_snv_data(self):
        prepared_muse_snv = prepare_muse_snv(
            get_tcga("muse_snv", self.cancer)
        )
        prepared_muse_snv_data = SomaticMutation(prepared_muse_snv)
        return prepared_muse_snv_data

    @property
    def mutect2_snv_data(self):
        prepared_mutect2_snv = prepare_mutect2_snv(
            get_tcga("mutect2_snv", self.cancer)
        )
        prepared_mutect2_snv_data = SomaticMutation(prepared_mutect2_snv)
        return prepared_mutect2_snv_data

    @property
    def somaticsniper_snv_data(self):
        prepared_somaticsniper_snv = prepare_somaticsniper_snv(
            get_tcga("somaticsniper_snv", self.cancer)
        )
        prepared_somaticsniper_snv_data = SomaticMutation(prepared_somaticsniper_snv)
        return prepared_somaticsniper_snv_data

###################
## MUTATION DATA ##
###################

def prepare_muse_snv(raw_muse_snv_data):
    return raw_muse_snv_data

def prepare_mutect2_snv(raw_mutect2_snv_data):
    mutect2_snv_data = (raw_mutect2_snv_data
                        .rename(columns=dict(Sample_ID="sample_id",
                                             gene="gene_name",
                                             chrom="Chromosome",
                                             start="Start",
                                             end="End"))
                        .drop(["ref", "alt", "Amino_Acid_Change"], axis=1))
    return mutect2_snv_data

def prepare_somaticsniper_snv(raw_somaticsniper_snv_data):
    return raw_somaticsniper_snv_data

######################
## COPY NUMBER DATA ##
######################

def prepare_cnv(raw_cnv_data):
    cnv_data = (raw_cnv_data
                .rename(columns=dict(Chrom="Chromosome",
                                     sample="sample_id")))
    cnv_data["Chromosome"] = cnv_data["Chromosome"].apply(lambda x: "chr" + str(x))
    return cnv_data

def prepare_gistic(raw_gistic_data):
    # TODO: Fails with pancancer data
    gistic_data = raw_gistic_data.set_index("Gene Symbol").transpose()
    gistic_data.columns = map(lambda x: x.split(".")[0], gistic_data.columns)
    return gistic_data

def prepare_masked_cnv(raw_masked_cnv_data):
    masked_cnv_data = (raw_masked_cnv_data
                       .rename(columns=dict(Chrom="Chromosome",
                                            sample="sample_id")))
    masked_cnv_data["Chromosome"] = masked_cnv_data["Chromosome"].apply(lambda x: "chr" + str(x))
    return masked_cnv_data

#####################
## EXPRESSION DATA ##
#####################

def prepare_htseq_counts(raw_htseq_counts):
    htseq_counts = raw_htseq_counts.set_index("Ensembl_ID").transpose()
    htseq_counts.columns = map(lambda x: x.split(".")[0], htseq_counts.columns)
    return htseq_counts

def prepare_htseq_fpkm_uq(raw_htseq_fpkm_uq):
    htseq_fpkm_uq = raw_htseq_fpkm_uq.set_index("Ensembl_ID").transpose()
    htseq_fpkm_uq.columns = map(lambda x: x.split(".")[0], htseq_fpkm_uq.columns)
    return htseq_fpkm_uq
