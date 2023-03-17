from .generics import simulate_matrix


def simple_mutation(
    n_samples,
    n_vars,
    *args,
    dna_vaf_q=0.9,
    random_state=None,
    **kwargs
):
    mutation = simulate_matrix(n_samples, n_vars, "snv_", random_state, rv="uniform")
    mutation = mutation.rename_axis("sample_id", axis=0).rename_axis("snv_id", axis=1)
    # Binarize
    mutation = (mutation
                .where(lambda df: df > dna_vaf_q, other=0)
                .where(lambda df: df <= dna_vaf_q, other=1))
    
    return mutation

mutation_methods = {
    "simple_mutation": simple_mutation,
}
