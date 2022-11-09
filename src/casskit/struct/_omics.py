from abc import ABC

from ..struct._base import AssayData


class Assay(ABC):
    pass

class CopynumberMixin(Assay):
    pass

class ExpressionMixin(Assay):
    pass

class MutationMixin(Assay):
    pass

class Copynumber(CopynumberMixin, AssayData):
    pass

class Expression(ExpressionMixin, AssayData):
    pass

class Mutation(MutationMixin, AssayData):
    pass
