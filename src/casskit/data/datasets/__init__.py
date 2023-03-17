from pkgutil import extend_path

from .tcga_dataset import TCGADataSet


__all__ = ["TCGADataSet"]
__path__ = extend_path(__path__, __name__)
