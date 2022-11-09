import logging
import os

os.environ["CASSKIT_CACHE_DIR"] = "/scratch/users/tsilvers/ceqtl_selection"

from ._tiledb import make_data_db

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# https://xenabrowser.net/datapages/

__all__ = [
    "make_data_db",
]