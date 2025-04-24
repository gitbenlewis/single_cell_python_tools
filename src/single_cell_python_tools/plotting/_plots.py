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

# set up logging within the module (not the root logger)
import logging
__name__leaf = __name__.split('.')[-1]
logger = logging.getLogger("sctl.pl." + __name__leaf)


sc.settings.set_figure_params(dpi=80, facecolor='white')



def silhouette_score_n_plot(adata,leiden_res='unk',**parameters):
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

def silhouette_score_of_obs_key_n_plot(adata,obs_key='leiden',**parameters):
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

def plot_batch_obs_key_of_obs_key2_old(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count"):
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

def plot_batch_obs_key_of_obs_key2(adata,savetable=False,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count"):
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

def plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4)):
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