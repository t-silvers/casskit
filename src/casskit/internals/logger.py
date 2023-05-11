import os
from pathlib import Path
import logging


class CassKitLogManager:
    def __init__(self, log_dir=None):
        if log_dir is not None:
            self.log_dir = Path(log_dir)
        else:
            self.log_dir = Path(os.getenv('PKG_LOG', Path.home() / '.demo_pkg_logs'))

        self.log_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(filename=str(self.log_dir / 'demo_pkg.log'), level=logging.INFO)
        self._logger = logging.getLogger(__name__)

    def get_log_dir(self):
        return str(self.log_dir)

    def set_log_dir(self, new_dir):
        self.log_dir = Path(new_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        for handler in self._logger.handlers[:]:
            self._logger.removeHandler(handler)
        logging.basicConfig(filename=str(self.log_dir / 'demo_pkg.log'), level=logging.INFO)

    def get_logger(self):
        return self._logger
