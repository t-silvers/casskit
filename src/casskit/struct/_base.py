from abc import ABC, abstractmethod


class DataMixin:
    pass

class SparseData(DataMixin, ABC):
    pass

class DenseData(DataMixin, ABC):
    pass

class FeatureData(SparseData):
    pass

class SampleData(SparseData):
    pass

class AssayData(DenseData):
    pass
