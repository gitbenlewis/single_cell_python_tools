'''
the plotting subpackage
'''
from __future__ import annotations
from ._plots import * 
from . import _plots
__all__: list[str] = _plots.__all__.copy()