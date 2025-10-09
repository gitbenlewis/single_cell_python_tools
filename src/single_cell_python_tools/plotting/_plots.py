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



def silhouette_score_n_plot(
    adata: anndata.AnnData | None = None,
    leiden_res: str | float | None = 'unk',
    **parameters: Any
):
    '''
    adsctl.pl.silhouette_score_n_plot(adata,parameters,leiden_res='unk'):
    > assumes ledien clusteirng to subset cells
    > uses X_pca for silhoutte_scores
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden']
  
    '''
    ##################### sillhouette scoreing
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden'])
    adata.obs['silhoutte']=samples_silhoutte_scores.tolist()
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    logger.info(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    #print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of add this to adata.uns later')

    ##################### sillhouette scoreing  #####END

    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting 
    fig_PP2C_cluster_scores, (ax_final,UMAP_final,UMAP_sil,pca_leiden) = plt.subplots(nrows=1, ncols=4, figsize=(20,5), gridspec_kw={'wspace':0.4})

    UMAP_final=sc.pl.umap(adata, color='leiden',title=str('leiden')+' Avg.sil.='+str(silhouette_score_adata), ax=UMAP_final,#palette=sc.pl.palettes.vega_20_scanpy,
                          show=False)

    cluster_silhouette_score_list=cluster_silhouette_score_list=adata.obs.groupby('leiden')[f'silhoutte'].mean()


    #  cluster scores
    pre_scores=cluster_silhouette_score_list
    pre_y_pos = cluster_silhouette_score_list.index.tolist()
    ax_final.barh(pre_y_pos,pre_scores)
    #ax_final.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_final.set_yticks(pre_y_pos)
    ax_final.set_yticklabels(pre_y_pos) #new
    ax_final.invert_yaxis()  # labels read top-to-bottom
    ax_final.set_title('Cluster Silhoutte Scores')

    UMAP_sil=sc.pl.umap(adata, color=['silhoutte'], ax=UMAP_sil, show=False,#palette=sc.pl.palettes.vega_20_scanpy
                       )

    pca_leiden=sc.pl.pca(adata, color='leiden', ax=pca_leiden, show=False)

    #fig_PP2C_cluster_scores.savefig(dataset_figures_output_directory+'silscore.pdf')
    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting  #####END
    #return adata
    return 

def silhouette_score_of_obs_key_n_plot(
    adata: anndata.AnnData | None = None,
    obs_key: str | None = 'leiden',
    **parameters: Any
):
    '''
    adsctl.pl.silhouette_score_n_plot(adata,obs_key='leiden',**parameters):
    > assumes ledien clusteirng to subset cells
    > uses X_pca for silhoutte_scores
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden']
    obs_key='leiden'
  
    '''
    ##################### sillhouette scoreing
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs[obs_key])
    adata.obs[f'silhoutte_{obs_key}']=samples_silhoutte_scores.tolist()
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs[obs_key],)
    cluster_number=len(set(adata.obs[obs_key].tolist()))
    logger.info(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} number of {obs_key} groups') 
    #print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of add this to adata.uns later')

    ##################### sillhouette scoreing  #####END

    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting 
    fig_PP2C_cluster_scores, (ax_final,UMAP_final,UMAP_sil,pca_leiden) = plt.subplots(nrows=1, ncols=4, figsize=(20,5), gridspec_kw={'wspace':0.4})

    UMAP_final=sc.pl.umap(adata, color=obs_key,title=str(obs_key)+' Avg.sil.='+str(silhouette_score_adata), ax=UMAP_final,#palette=sc.pl.palettes.vega_20_scanpy,
                          show=False)

    cluster_silhouette_score_list=cluster_silhouette_score_list=adata.obs.groupby(obs_key)[f'silhoutte_{obs_key}'].mean()


    #  cluster scores
    pre_scores=cluster_silhouette_score_list
    pre_y_pos = cluster_silhouette_score_list.index.tolist()
    ax_final.barh(pre_y_pos,pre_scores)
    #ax_final.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_final.set_yticks(pre_y_pos)
    ax_final.set_yticklabels(pre_y_pos) #new
    ax_final.invert_yaxis()  # labels read top-to-bottom
    ax_final.set_title('Cluster Silhoutte Scores')

    UMAP_sil=sc.pl.umap(adata, color=[f'silhoutte_{obs_key}'], ax=UMAP_sil, show=False,#palette=sc.pl.palettes.vega_20_scanpy
                       )

    pca_leiden=sc.pl.pca(adata, color=obs_key, ax=pca_leiden, show=False)

    #fig_PP2C_cluster_scores.savefig(dataset_figures_output_directory+'silscore.pdf')
    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting  #####END
    #return adata
    return


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

def plot_batch_obs_key_of_obs_key2(
    adata: anndata.AnnData | None = None,
    batch_obs_key: str = 'batch',
    obs_key2: str = "leiden",
    flavor: str = "pct_count",
    figsize: tuple = (10, 4),
    savetable: bool = False,
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
    import os
    if savefig:
        os.makedirs(output_dir+output_prefix+'/figures/', exist_ok=True)
    if savetable:
        os.makedirs(output_dir+output_prefix+'/tables/', exist_ok=True)
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
        if savetable==True:
            df_norm.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_norm.csv')
            df.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_counts.csv')
    if flavor=="pct":
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax1=df_norm.plot.bar(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
        if savetable==True:
            df_norm.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_norm.csv')
    if flavor=="count":
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax2=df.plot.bar(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
        if savetable==True:
            df.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_counts.csv')
    if flavor=="pct_count_barh":
        df_norm_sorted = df_norm.sort_index( ascending=False)
        df_sorted = df.sort_index(ascending=False)
        fig1, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize,sharey=True)
        ax1=df_norm_sorted.plot.barh(stacked=True,ax=axes[0]).legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),  borderaxespad=0)
        # add x-axis label to the left chart
        axes[0].set_xlabel("Fraction of cells")
        ax2=df_sorted.plot.barh(stacked=True,ax=axes[1]).legend().set_visible(False)
        # add x-axis label to the right chart
        axes[1].set_xlabel("Number of cells")
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
        if savetable==True:
            df_norm_sorted.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_norm.csv')
            df_sorted.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_counts.csv')
    if flavor=="pct_barh":
        df_norm_sorted = df_norm.sort_index( ascending=False)
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax1=df_norm_sorted.plot.barh(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
        if savetable==True:
            df_norm_sorted.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_norm.csv')
    if flavor=="count_barh":
        df_sorted = df.sort_index(ascending=False)
        fig1, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize)
        ax2=df_sorted.plot.barh(stacked=True,ax=axes).legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        if savefig==True:
            fig1.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'.pdf')
        if savetable==True:
            df_sorted.to_csv(output_dir+output_prefix+'/tables/'+output_prefix+'crosstab_'+obs_key2+'_'+batch_obs_key+'_counts.csv')
    # make sure to show the plot
    plt.show()
    return df, df_norm

####################################################################################################################

####################################################################################################################


#### plots for adata distributions #################################################################################################################

def plot_adata_raw_and_X(
    adata: anndata.AnnData | None = None,
    max_value_mask: float| None = 3e4,
    savefig: bool| None = False,
    output_dir: str| None = './project/',
    output_prefix: str| None = "dataset_",
    **kwargs
):
    '''plot adata raw counts and adata.X counts side by side'''
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(1, 2, figsize=(10, 6))
    # set figure title
    fig.suptitle('Distribution of Feature Values in layers of adata object', fontsize=16)
    row_totals_raw = adata.raw.X.sum(axis=1).A1
    row_totals_X = adata.X.sum(axis=1).A1
    ## plot adata.raw
    raw_mask = row_totals_raw < max_value_mask
    ax[0].hist(row_totals_raw[raw_mask], bins=100, color='blue', alpha=0.7)
    stats_string_0=f'min: {row_totals_raw.min():.1f}\nmax: {row_totals_raw.max():.1f}\nmean: {row_totals_raw.mean():.1f}\nstd: {row_totals_raw.std():.1f}\n25%: {np.percentile(row_totals_raw, 25):.1f}\n50%: {np.percentile(row_totals_raw, 50):.1f}\n75%: {np.percentile(row_totals_raw, 75):.1f}\n99.9%: {np.percentile(row_totals_raw, 99.9):.1f}\ncount: {row_totals_raw.shape[0]:.0f}'
    ax[0].set_title(f'adata.raw feature matrix')
    ax[0].set_xlabel(f'Feature values per observation\n{stats_string_0}')
    ax[0].set_ylabel('Number of observations')
    ## plot adata.X
    X_mask = row_totals_X < max_value_mask
    ax[1].hist(row_totals_X[X_mask], bins=100, color='green', alpha=0.7)
    stats_string_1=f'min: {row_totals_X.min():.1f}\nmax: {row_totals_X.max():.1f}\nmean: {row_totals_X.mean():.1f}\nstd: {row_totals_X.std():.1f}\n25%: {np.percentile(row_totals_X, 25):.1f}\n50%: {np.percentile(row_totals_X, 50):.1f}\n75%: {np.percentile(row_totals_X, 75):.1f}\n99.9%: {np.percentile(row_totals_X, 99.9):.1f}\ncount: {row_totals_X.shape[0]:.0f}'
    ax[1].set_title(f'adata.X feature matrix')
    ax[1].set_xlabel(f'Feature values per observation\n{stats_string_1}')
    ax[1].set_ylabel('Number of observations')
    plt.tight_layout()
    plt.show()
    if savefig==True:
        os.makedirs(output_dir+output_prefix+'/figures/', exist_ok=True)
        fig.savefig(output_dir+output_prefix+'/figures/'+output_prefix+'adata_raw_and_X_dist.pdf')    
    return fig, ax


##### END plots for adata distributions ###############################################################################################################




# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]