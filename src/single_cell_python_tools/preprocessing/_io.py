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

# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]