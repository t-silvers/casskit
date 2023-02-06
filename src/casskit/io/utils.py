from functools import wraps
from pathlib import Path
import pkg_resources
import platform
import subprocess
from typing import Any, Callable, Dict, List, Union
import warnings

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

            if (Path(cache).exists() and Path(cache).stat().st_size > 0):
                print(f"Loading from cache: {cache}")
                data = self.read_cache(cache)
                
            else:
                print(f"Caching to disk: {cache}")
                Path(cache).parent.mkdir(exist_ok=True, parents=True)
                data = f(self, *args, **kwargs)
                if (data is not None and data.empty is False):
                    self.write_cache(data, cache)
                else:
                    warnings.warn(f"Data is empty. Not writing to cache: {cache}")

        else:
            data = f(self, *args, **kwargs)

        return data
    return wrapper

def check_package_version(package: str, version: str = None) -> bool:
    """Check if package is installed and at least a certain version."""
    try:
        if package == "python":
            if version is None:
                return True
            else:
                return pkg_resources.parse_version(
                    platform.python_version()
                ) >= pkg_resources.parse_version(version)
        
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

def translate_pcawg_cancer(s):
    """Harmonize project cancer types"""
    PCAWG_CODES = {
        'BLCA-US': 'TCGA-BLCA', 'BRCA-US': "TCGA-BRCA", 'OV-AU': "TCGA-OV",
        'PAEN-AU': "TCGA-PDAC", 'PRAD-CA': "TCGA-PRAD", 'PRAD-US': "TCGA-PRAD",
        'RECA-EU': "TCGA-KIRC|TCGA-KIRP", 'SKCM-US': "TCGA-SKCM",
        'STAD-US': "TCGA-STAD", 'THCA-US': "TCGA-THCA", 'KIRP-US': "TCGA-KIRP",
        'LIHC-US': "TCGA-LIHC", 'PRAD-UK': "TCGA-PRAD", 'LIRI-JP': "TCGA-LIHC",
        'PBCA-DE': np.nan, 'CESC-US': "TCGA-CESC", 'PACA-AU': "TCGA-PDAC",
        'PACA-CA': "TCGA-PDAC", 'LAML-KR': "TCGA-LAML", 'COAD-US': "TCGA-COAD",
        'ESAD-UK': "TCGA-ESCA", 'LINC-JP': "TCGA-LIHC", 'LICA-FR': "TCGA-LIHC",
        'CLLE-ES': np.nan, 'HNSC-US': "TCGA-HNSC", 'EOPC-DE': "TCGA-PRAD",
        'BRCA-UK': "TCGA-BRCA", 'BOCA-UK': np.nan, 'MALY-DE': "TCGA-DLBC",
        'CMDI-UK': np.nan, 'BRCA-EU': "TCGA-BRCA", 'ORCA-IN': np.nan,
        'BTCA-SG': "TCGA-CHOL", 'SARC-US': "TCGA-SARC", 'KICH-US': "TCGA-KICH",
        'MELA-AU': "TCGA-SKCM", 'DLBC-US': "TCGA-DLBC", 'GACA-CN': "TCGA-STAD",
        'PAEN-IT': "TCGA-PDAC", 'GBM-US': "TCGA-GBM", 'KIRC-US': "TCGA-KIRC",
        'LAML-US': "TCGA-LAML", 'LGG-US': "TCGA-LGG", 'LUAD-US': "TCGA-LUAD",
        'LUSC-US': "TCGA-LUSC", 'OV-US': "TCGA-OV", 'READ-US': "TCGA-READ",
        'UCEC-US': "TCGA-UCEC",
    }
    return s.rename("pcawg_id").replace(PCAWG_CODES).str.split("|")