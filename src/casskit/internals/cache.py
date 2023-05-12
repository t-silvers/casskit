import os
from pathlib import Path
import shutil

try:
    from pypath.share import settings
    can_pypath = True
except ImportError:
    can_pypath = False


class CassKitCacheManager:
    def __init__(self, cache_dir=None):
        if cache_dir is not None:
            self.cache_dir = Path(cache_dir)
        else:
            self.cache_dir = Path(os.getenv('PKG_CACHE', Path.home() / '.demo_pkg_casche'))

        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get_cache_dir(self):
        return str(self.cache_dir)

    def set_cache_dir(self, new_dir):
        self.cache_dir = Path(new_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # If using pypath, match its cache directory
        # TODO: This should only be implemented if the use is using pypath via cass-kit.
        #       The current implementation affects already installed pypath.
        if can_pypath:
            settings.setup(cachedir=Path(new_dir) / 'pypath',
                        pickle_dir=Path(new_dir) / 'pypath' / 'pickles')

    def clear_cache(self):
        shutil.rmtree(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def cache_size(self):
        return sum(f.stat().st_size for f in self.cache_dir.glob('**/*') if f.is_file())