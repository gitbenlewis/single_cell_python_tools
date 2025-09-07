# module level import libraries
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
import anndata
from typing import Any, Dict, Optional, List

# set up logging within the module (not the root logger)
import logging
__name__leaf = __name__.split('.')[-1]
logger = logging.getLogger("sctl.pl." + __name__leaf)


sc.settings.set_figure_params(dpi=80, facecolor='white')




####################################################################################################################

####################################################################################################################

def plot_percent_obs_key2_per_batch_obs_key(
    adata: anndata.AnnData | None = None,
    batch_obs_key='batch',
    obs_key2="leiden",
    figsize=(10,4),
    savefig=False,
    output_dir='./project/',
    output_prefix="dataset_",
):
    '''
     df_norm=adsctl.pl.plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4)
     This produce one column of individual bar charts (one chart for each catagory in batch_obs_key='batch'") # batch_obs_key="sample_ID is good to use
     each bar chart show percentage of cells in "batch" assigned to obs_key2="leiden"
    returns df_col_norm
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    obs_key2_Xtab_batch_obs_key_df=pd.crosstab(adata.obs[obs_key2],adata.obs[batch_obs_key] )

    df=obs_key2_Xtab_batch_obs_key_df

    df_col_norm=pd.DataFrame()
    for i in df.columns:
        df_col_norm[i]=list(map(lambda x:x/df[i].sum(axis=0),df[i]))
    #print(df_col_norm)
    fig1, axes = plt.subplots(nrows=df_col_norm.shape[1], ncols=1,figsize=figsize)
    ax_n=0
    #obs_key2_groups=np.arange(len(df_col_norm.index))
    obs_key2_groups=adata.obs[obs_key2].cat.categories.tolist()
    for i in df_col_norm.columns:
        #ax=df_col_norm[df_col_norm.columns[ax_n]].plot.bar(ax=axes[ax_n]).legend().set_visible(True)
        axes[ax_n].barh(obs_key2_groups,df_col_norm[df_col_norm.columns[ax_n]].tolist())
        axes[ax_n].set_title(f' {df_col_norm.columns[ax_n]}')
        axes[ax_n].set_yticks(obs_key2_groups, labels=obs_key2_groups)
        axes[ax_n].invert_yaxis() 

        for bars in axes[ax_n].containers:
            axes[ax_n].bar_label(bars, label_type='center',  fmt='%.2g',padding=30,)
        ax_n=ax_n+1
    if savefig==True:
        fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    return df_col_norm 

####################################################################################################################




# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]