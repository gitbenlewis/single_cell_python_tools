## module imports
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import anndata
from typing import Any, Dict, Optional, List


# set up logging within the module (not the root logger)
import logging
__name__leaf = __name__.split('.')[-1]
logger = logging.getLogger("sctl.pp." + __name__leaf)


def reset_cellxgene_var_names(adata: anndata.AnnData,
                              feature_name_col: str = 'feature_name',
                            **parameters: Any) -> None:
    """
    Reset var names for cellxgene datasets.

    This function modifies the AnnData object in place. It creates a new column
    in `adata.var` called 'feature_name_unique' that appends the value in the
    index column (default 'gene_id') if the value in the name column (default
    'gene_name') is not unique. The var names of the AnnData object are then set
    to this new column.
    index column (default 'gene_id') if the value in the name column (default
    'gene_name') is not unique. The var names of the AnnData object are then set
    to this new column.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to modify.
    index_col : str, optional
        The column in `adata.var` to use as the unique identifier (default is 'gene_id').
    name_col : str, optional
        The column in `adata.var` to use as the feature name (default is 'gene_name').

    Returns
    -------
    None
        The function modifies the AnnData object in place and does not return anything.
    """
    if feature_name_col not in adata.var.columns:
        raise ValueError(f"'{feature_name_col}'  must be columns in adata.var")

    # Create a new column 'feature_name_unique'
    adata.var[f'{feature_name_col}_unique'] = [
        f"{name}_{idx}" if duplicated else name
        for name, idx, duplicated in zip(
            adata.var[feature_name_col],
            adata.var.index,
            adata.var[feature_name_col].duplicated(keep=False)
        )
    ]

    # Set the new var names
    adata.var.reset_index(inplace=True)
    adata.var.set_index(f'{feature_name_col}_unique', inplace=True)

    # Update raw if it exists
    if adata.raw is not None:
        adata.raw = adata.copy()  # reset raw to adata with new var_names

    logger.info(f"Var_names reset using column '{feature_name_col}' \nwith uniqueness ensured by appending index where necessary.")


import numpy as np
import scipy.sparse as sp
import anndata as ad
from collections.abc import Sequence
from typing import Any

def _as_float32(matrix):
    if matrix is None:
        return matrix
    if sp.issparse(matrix):
        return matrix.astype(np.float32)
    return np.asarray(matrix, dtype=np.float32)

def downcast_layers_to_float32(
    adata: ad.AnnData,
    *,
    include_X: bool = True,
    include_raw: bool = False,
    exclude_layers: Sequence[str] | None = None,
    **kwargs: Any,
) -> None:
    """Convert AnnData layers (and optionally X/raw) to float32 in place."""
    _ = kwargs  # intentionally unused; allows extra config to pass through
    skip = set(exclude_layers or ())
    for layer_key in list(adata.layers.keys()):
        if layer_key in skip:
            continue
        adata.layers[layer_key] = _as_float32(adata.layers[layer_key])
    if include_X:
        adata.X = _as_float32(adata.X)
    if include_raw and adata.raw is not None:
        raw = adata.raw.to_adata()
        raw.X = _as_float32(raw.X)
        adata.raw = raw

# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]