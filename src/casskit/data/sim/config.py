DEFAULT_CFG = {
    "order": ["copynumber", "mutation", "grn", "expression"],
    "copynumber": {
        "method": "chromswap_augmented",
        "kwargs": {},
        "alias": "cnv_data",
    },
    "mutation": {
        "method": "simple_mutation",
        "kwargs": {},
        "alias": "mutect2_snv_data",
    },
    "expression": {
        "method": "from_grn",
        "kwargs": {},
        "alias": "htseq_counts_data",
    },
    "grn": {
        "kwargs": {},
    },
}

BUILDER_ALIASES = {
    "copynumber": "CopyNumberVariationBuilder",
    "mutation": "SomaticMutationBuilder",
    "expression": "MessengerRNABuilder",
}