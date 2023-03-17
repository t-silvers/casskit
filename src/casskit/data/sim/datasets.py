import copy
import sys
sys.path.append("/Users/thomassilvers/GitHub/")
from casskit_data.multiomic import MultiOmic

from . import builders
from .builders.base import OmicSimulator
from .builders.grn import GRNgineer
from .config import DEFAULT_CFG, BUILDER_ALIASES


def simulate_tcga(
    n_samples,
    n_vars,
    *args,
    config=DEFAULT_CFG,
    copynumber_method=None,
    mutation_method=None,
    expression_method=None,
    **kwargs,
):
    # Initiate simulated data set class and builder class
    sim_multiomic = MultiOmic()
    omic_simulator = OmicSimulator(n_samples, n_vars)
    
    # Parse config and build
    config = copy.copy(config)
    if copynumber_method is not None:
        config["copynumber"]["method"] = copynumber_method
    if mutation_method is not None:
        config["mutation"]["method"] = copynumber_method
    if expression_method is not None:
        config["expression"]["method"] = copynumber_method
    
    kwargs_rollover = kwargs
    builder_order = config.get("order", DEFAULT_CFG["order"])
    design_builders = dict()
    for builder_name in builder_order:
        builder_cfg = config.get(builder_name, {})
        builder_kwargs = builder_cfg.get("kwargs", {}).copy()
        builder_kwargs.update(kwargs_rollover)
        
        if builder_name == "grn":
            grn, design = GRNgineer.from_datasets(design_builders, n_vars, **builder_kwargs)
            
            # Allow GRN output to rollover to expression builder
            kwargs_rollover.update(dict(grn=grn, design=design))
        
        else:
            builder_method = builder_cfg.get("method", None)
            if builder_method is None:
                raise ValueError(f"Method not specified for {builder_name}")
            builder_class = getattr(builders, BUILDER_ALIASES[builder_name])
            builder = builder_class(**builder_kwargs)
            omic_simulator.simulate_omics(builder, builder_method, **builder_kwargs)
            assert omic_simulator.method == builder_method
            builder_alias = builder_cfg.get("alias", builder_name)
            setattr(sim_multiomic, builder_alias, omic_simulator.simulated_omics.data)

            # Allow output to rollover to GRN builder
            design_builders.update({builder_name: omic_simulator.simulated_omics.data})

    return sim_multiomic, grn