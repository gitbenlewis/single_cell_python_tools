"""  
Bens single cell analysis tools ... mostly scanpy wrappers
(importable as ``import single_cell_python_tools as sctl``)
Exposes
-------
pp  - preprocessing sub-namespace
pl  - plotting sub-namespace
tl  - general-purpose tools
DATASET_class - high-level class

"""
# logging
from __future__ import annotations
import logging
import sys

# suppress future warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# ───────────────────────────────
# LOGGING  (one handler, one place)
# ───────────────────────────────
_LOGGER_NAME = "sctl"   
_LOG_LEVEL   = logging.INFO   # or DEBUG

def _configure_root_logger() -> None:
    logger = logging.getLogger(_LOGGER_NAME)
    logger.setLevel(_LOG_LEVEL)

    # Add exactly one handler if none present (import-safe)
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        handler = logging.StreamHandler(stream=sys.stdout)
        fmt = "%(asctime)s | %(name)s | %(levelname)-s | %(message)s"
        handler.setFormatter(logging.Formatter(fmt))
        logger.addHandler(handler)

    # Prevent the log records from being duplicated by the *global* root
    logger.propagate = False
# run at import time
_configure_root_logger() 

# ──────────────────────────────────────────────────────────────
# Public API – 
# ──────────────────────────────────────────────────────────────
from . import preprocessing as pp          # sub-namespaces
from . import plotting      as pl
from . import tools         as tl

from . import DATASET_class  # high-level class
from .DATASET_class import DATASET_class # to access the class as sctl.DATASET_class(...)

__all__: list[str] = [
    "pp",                   # sctl.pp
    "pl",                   # sctl.pl
    "tl",                   # sctl.tl
    "DATASET_class",   # sctl.DATASET_class
]