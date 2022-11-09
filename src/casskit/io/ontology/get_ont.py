from ...io.ontology._cosmic import get_cosmic


def get_ont(ontology_id):
    """Get ontology data from local cache."""

    match ontology_id:
        
        case "cosmic":
            return get_cosmic()
                
        case _:
            raise ValueError(f"Feature {ontology_id} not found.")