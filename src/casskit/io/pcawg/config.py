from collections import namedtuple


PCAWGData = namedtuple("pcawg", ["stem", "omic"])

PCAWG_XENA_DATASETS = {
    "rnaseq": PCAWGData("tophat_star_fpkm_uq.v2_aliquot_gl.sp.log", "rnaseq"),
    "copynumber": PCAWGData("20170119_final_consensus_copynumber_sp", "copynumber"),
    "phenotype": PCAWGData("project_code_sp", "phenotype"),
}

