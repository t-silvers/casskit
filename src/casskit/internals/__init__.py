from .cache import CassKitCacheManager
from .logger import CassKitLogManager


# Set default internal managers for cass-kit
_cache_manager = CassKitCacheManager()
_logger_manager = CassKitLogManager()
_logger = _logger_manager.get_logger()


def set_cache_dir(new_dir):
    _cache_manager.set_cache_dir(new_dir)

def get_cache_dir():
    return _cache_manager.get_cache_dir()

def clear_cache():
    _cache_manager.clear_cache()

def cache_size():
    return _cache_manager.cache_size()

def set_log_dir(new_dir):
    _logger_manager.set_log_dir(new_dir)

def get_log_dir():
    return _logger_manager.get_log_dir()

def get_logger():
    return _logger