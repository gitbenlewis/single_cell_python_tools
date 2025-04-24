"""  Bens single cell analysis tools ... mostly scanpy wrappers"""
# logging
from __future__ import annotations
import logging
import sys

_LOGGER_NAME = "sctl"   
_LOG_LEVEL   = logging.INFO   # or DEBUG

def _configure_root_logger() -> None:
    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(_LOG_LEVEL)

    # Add exactly one handler if none present (import-safe)
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        handler = logging.StreamHandler(stream=sys.stdout)
        fmt = "%(asctime)s | %(name)s | %(levelname)-8s | %(message)s"
        handler.setFormatter(logging.Formatter(fmt))
        logger.addHandler(handler)

    # Prevent the log records from being duplicated by the *global* root
    logger.propagate = False

_configure_root_logger()           # run at import time

from . import preprocessing as pp
from . import plotting as pl
from . import tools as tl

# from . import preprocessing *
# from . import plotting *
# from . import tools *
