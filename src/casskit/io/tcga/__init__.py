from . import get as get_tcga

__path__ = __import__('pkgutil').extend_path(__path__, __name__)
__all__ = ["get_tcga"]