import logging
import pytest
import re


def setup_logging(level):
    if level is None:
        level = logging.WARN
    elif re.match(r"\d+", level):
        level = int(level)
    logging.basicConfig(level=level)

@pytest.fixture
def datadir(tmpdir_factory):
    """Setup base cache directory for a test"""
    p = tmpdir_factory.mktemp(".cache")
    return p
