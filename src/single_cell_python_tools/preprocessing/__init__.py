"""
the Pre-processing helpers sub package 
(imported as ``sctl.pp``).
Sub-modules:
    _qc              – quality-control helpers
    _transform_data  – normalisation / HVG / PCA
    _clustering      – Leiden & silhouette utilities
"""
from __future__ import annotations

# Import the three implementation files as public sub-namespaces
from . import _qc
from . import _transform_data
from . import _clustering

# Re-export their user-facing symbols so that
# `sctl.pp.basic_filtering()` works.
from ._qc              import *  
from ._transform_data  import * 
from ._clustering      import *  

__all__: list[str] = []    
for _m in (_qc, _transform_data, _clustering):
    __all__.extend(getattr(_m, "__all__", []))