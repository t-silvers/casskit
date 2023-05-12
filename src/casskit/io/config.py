import os
from pathlib import Path


DEFAULT_CACHE = Path("~/.cache").expanduser()

try:
    user_cache = os.environ["CASSKIT_CACHE_DIR"]
    CACHE_DIR = Path(user_cache)
except KeyError:
    CACHE_DIR = DEFAULT_CACHE

HEADER = {
    "User-Agent": ("Mozilla/5.0 (Macintosh;"
                   "Intel Mac OS X 10_14_6)"
                   "AppleWebKit/605.1.15 (KHTML, like Gecko) "
                   "Version/14.1.2 Safari/605.1.15"),
    "Accept": ("text/html,"
               "application/xhtml+xml,"
               "application/xml;q=0.9,"
               "image/webp,"
               "*/*;q=0.8")
}