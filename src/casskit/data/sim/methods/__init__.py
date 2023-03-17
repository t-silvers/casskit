from pkgutil import extend_path

from .copynumber import copynumber_methods
from .expression import expression_methods
from .mutation import mutation_methods


simulation_methods = dict(
    CopyNumberVariation=copynumber_methods,
    MessengerRNA=expression_methods,
    SomaticMutation=mutation_methods,
)

__all__ = ["simulation_methods"]
__path__ = extend_path(__path__, __name__)
