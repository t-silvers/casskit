import logging
import os

os.environ["CASSKIT_CACHE_DIR"] = "/scratch/users/tsilvers/ceqtl_selection"

# from ._tiledb import DataDB

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# __all__ = [
#     "DataDB",
# ]