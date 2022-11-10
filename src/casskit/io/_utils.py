from functools import wraps
from pathlib import Path
import pkg_resources
import subprocess
from typing import Any, Callable, Dict, List, Union

import pandas as pd


def cache_on_disk(f: Callable) -> Callable:
    """Cache function output on disk.
    
    Light-weight caching decorator that caches class function output on disk.
    This decorator will only work on functions that return a single pandas
    DataFrame object.    
    """
    @wraps(f)
    def wrapper(self, *args, **kwargs):
        if hasattr(self, "path_cache"):
            cache = self.path_cache

            if Path(cache).exists():
                print(f"Loading from cache: {cache}")
                data = self.read_cache(cache)
                
            else:
                print(f"Caching to disk: {cache}")
                Path(cache).parent.mkdir(exist_ok=True, parents=True)
                data = f(self, *args, **kwargs)
                self.write_cache(data, cache)

        else:
            data = f(self, *args, **kwargs)

        return data
    return wrapper

def check_package_version(package: str, version: str = None) -> bool:
    """Check if package is installed and at least a certain version."""
    try:
        pkg_distrib = pkg_resources.get_distribution(package)
        if version is not None:
            return pkg_resources.parse_version(pkg_distrib.version) >= pkg_resources.parse_version(version)
    
    except pkg_resources.DistributionNotFound:
        return False

def wrap_rcall(script: str,
              params: Dict[Union[str, Path], pd.DataFrame],
              output: List[Path],
              check: bool = True) -> None:
    """Wrap R call to script."""
    def _args():
        for arg, param in params.items():
            param_ = param
            if isinstance(param, Path):
                param_ = param.as_posix()
            yield f'--{arg} {param_}'
    
    cmd = f"Rscript {script} {' '.join(_args())}"
    subprocess.run(cmd, shell=True, check=check)
