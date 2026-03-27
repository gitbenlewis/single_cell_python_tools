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
def plot_batch_obs_key_of_obs_key2_old(
    adata: anndata.AnnData | None = None,
    batch_obs_key: str = 'batch',
    obs_key2: str = "leiden",
    flavor: str = "pct_count",
    figsize: tuple[int, int] = (10, 4),
    savefig: bool = False,
    output_dir: str = './project/',
    output_prefix: str = "dataset_"
):
    '''
    df, df_norm=adsctl.pl.plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count")
    makes two side by side bar charts, each bar is a batch_obs_key='batch category and each bar is stacked and colored by obs_key2="leiden"
    left bar chart is fraction on y -axis 
    right  bar chart is obs/cell count on y -axis 
    flavor="pct_count"  >>> both charts
    flavor="pct"  >>> only pct chart
    flavor="count"  >>> only count chart

    flavor="pct_count_barh"  >>> both charts
    flavor="pct_barh"  >>> only pct chart
    flavor="count_barh"  >>> only count chart
    returns df, df_col_norm
    '''
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    df_norm=pd.crosstab(adata.obs[obs_key2],adata.obs[batch_obs_key], normalize='index')
    #print(df_norm)
    df=pd.crosstab(adata.obs[obs_key2],adata.obs[batch_obs_key] )
    #print(df)
    if flavor=="pct_count":
        fig1, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize)
        ax1=df_norm.plot.bar(stacked=True,ax=axes[0]).legend().set_visible(False)
        ax2=df.plot.bar(stacked=True,ax=axes[1]).legend(loc='upper right')
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="pct":
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax1=df_norm.plot.bar(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="count":
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax2=df.plot.bar(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="pct_count_barh":
        df_norm_sorted = df_norm.sort_index( ascending=False)
        df_sorted = df.sort_index(ascending=False)
        fig1, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize)
        ax1=df_norm_sorted.plot.barh(stacked=True,ax=axes[0]).legend().set_visible(False)
        ax2=df_sorted.plot.bar(stacked=True,ax=axes[1]).legend(loc='upper right')
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="pct_barh":
        df_norm_sorted = df_norm.sort_index( ascending=False)
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax1=df_norm_sorted.plot.barh(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    if flavor=="count_barh":
        df_sorted = df.sort_index(ascending=False)
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax2=df_sorted.plot.barh(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
    return df, df_norm

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