## module imports
import scanpy as sc
import pandas as pd
import anndata
from typing import Any, Dict, Optional, List

# set up logging within the module (not the root logger)
import logging
__name__leaf = __name__.split('.')[-1]
logger = logging.getLogger("sctl.pp." + __name__leaf)



def leiden_clustering(
    adata: anndata.AnnData | None = None,
    number_of_neighbors: int = 10,
    number_of_PC: int = 40,
    leiden_res: float = 1,
    key_added: str = 'leiden',
    **parameters: Any
):
    '''
    #### code
    leiden_clustering(adata, number_of_neighbors=10,number_of_PC=40, leiden_res=1,key_added='leiden', **parameters)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_res,key_added=key_added)
    sc.pl.umap(adata, color=[key_added] )
    '''
    sc.pp.neighbors(adata,n_neighbors=number_of_neighbors, n_pcs=number_of_PC)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_res,key_added=key_added)
    sc.pl.umap(adata, color=[key_added] )
    return


def rename_leiden_clusters_old(
    adata: anndata.AnnData | None = None,
    rename_cluster: bool = False,
    new_cluster_names: List[str] | None = None,
    new_obs_key: str = 'Cell_Clusters',
    **parameters: Any
):
    '''
    rename_leiden_clusters_old(adata,parameters=None,rename_cluster=False,new_cluster_names=None,new_obs_key='Cell_Clusters'):
    '''
    if rename_cluster==True:
        cluster_numbers_len=len(set(adata.obs['leiden'].tolist()))
        if cluster_numbers_len==len(new_cluster_names):
            adata.obs[new_obs_key]=adata.obs['leiden']
            adata.rename_categories(new_obs_key, new_cluster_names)
    return

def rename_leiden_clusters(
    adata: anndata.AnnData | None = None,
    rename_cluster: bool = False,
    new_cluster_names: List[str] | None = None,
    new_obs_key: str = 'Cell_Clusters_Named',
    reorder_cluster_names: bool = False,
    new_cluster_names_order: List[str] | None = None,
    **parameters: Any
):
    '''
    rename_leiden_clusters_old(adata,parameters=None,rename_cluster=False,new_cluster_names=None,new_obs_key='Cell_Clusters_Named'):
    '''
    if rename_cluster==True:
        leiden_number_list=adata.obs['leiden'].value_counts().index.tolist()
        #new_cluster_names=adata.uns["parameters"]['new_cluster_names']
        ## make a dictionary from leiden_number_list and new_cluster_names_list
        cluster_name_dict=dict(zip(leiden_number_list,new_cluster_names))
        adata.obs[new_obs_key]=adata.obs['leiden']
        adata.obs[new_obs_key] = adata.obs[new_obs_key].map(cluster_name_dict).fillna('Unknown').astype('category')
        if reorder_cluster_names:
            adata.obs[new_obs_key] = adata.obs[new_obs_key].cat.reorder_categories(new_cluster_names_order, ordered=True)
    return

def rename_leiden_clusters(
    adata: anndata.AnnData | None = None,
    rename_cluster: bool = False,
    new_cluster_names: List[str] | None = None,
    new_obs_key: str = 'Cell_Clusters_Named',
    reorder_cluster_names: bool = False,
    new_cluster_names_order: List[str] | None = None,
    **parameters: Any
):
    '''
    Rename Leiden clusters using a provided list of new cluster names.
    Parameters:
    - adata: AnnData object
    - rename_cluster: Boolean flag to perform renaming
    - new_cluster_names: List of cluster name strings in the order matching leiden cluster IDs
    - new_obs_key: Name of the new column to store cluster names
    - reorder_cluster_names: Whether to set the categories as an ordered Categorical
    - new_cluster_names_order: Explicit order for reordering categories
    '''
    if rename_cluster:
        leiden_number_list = adata.obs['leiden'].value_counts().index.tolist()
        cluster_name_dict = dict(zip(leiden_number_list, new_cluster_names))

        # Map and check for missing mappings
        mapped = adata.obs['leiden'].map(cluster_name_dict)

        if mapped.isna().any():
            categories = new_cluster_names + ['Unknown']
            mapped = mapped.fillna('Unknown')
        else:
            categories = new_cluster_names

        adata.obs[new_obs_key] = pd.Categorical(mapped, categories=categories, ordered=reorder_cluster_names)
        # reorder categories based on input list
        adata.obs[new_obs_key] = adata.obs[new_obs_key].cat.reorder_categories(categories, ordered=True)
        # reorder based on a specific order if provided
        if reorder_cluster_names and new_cluster_names_order is not None:
            adata.obs[new_obs_key] = adata.obs[new_obs_key].cat.reorder_categories(new_cluster_names_order, ordered=True)
    return


######################################################################### 
############ below were copied from old code with minimal name change

# needs a rewrtie looks ugly and is slow but works
    
# This file silhouette_of_leiden_clustering.py
#contains 3 funcitons 
#See doc strings for descriptions
# MD_leiden_cluster_sil_score(adata, leiden_res=1, n_jobs=8,savefig=False,output_dir="./figures/",output_prefix="adata"):
# MD_silhouette_walk_Largest_drop(adata,leiden_res_start=0.025,leiden_res_end=2,leiden_res_step=0.025,print_per_step=False, n_jobs=8,savetable=False,savefig=False,output_dir="./figures/",output_prefix="adata"):
# MD_silhouette_walk_4_Largest_drops(adata,leiden_res_start=0.025,leiden_res_end=2,leiden_res_step=0.025,print_per_step=False,n_jobs=8,savetable=False,savefig=False,output_dir="./figures/",output_prefix="adata"):
#
#########################################################################


""
def leiden_cluster_sil_score(
    adata: anndata.AnnData | None = None,
    leiden_res: float = 1.0,
    n_jobs: int = 1,
    savefig: bool = False,
    output_dir: str = "./figures/",
    output_prefix: str = "adata",
    **parameters: Any
):
    """
    defualt values :
    leiden_cluster_sil_score(
    leiden_res=1, n_jobs=8,savefig=False,output_dir="./figures/",output_prefix="adata")
    
    This function takes an adata object (required arguement)
    and performs leiden at a single leiden resolution value defualt value of 1 (leiden_res=1)
    
    It also calculates:
    Calculates average sillhoutte score for all observations at leiden_res
    Calculates average sillhoutte score of the observations in each cluster at leiden_res
    
    It then plots 4 plots side by side they are (left to right)
    1) Umap plot colored by the leiden clusters. 
    The plot title contains the leiden resolution value used and average sillhoutte score for all observations
    2) Bar graph of sillhouette score of each cluster 
    3) Umap plot colored by sillhouette score of each data point (cell) 
    4) PCA plot of first two PCAs colored by leiden clusters
    """
    #############################################################################   Import everything needed 
    import os
    import numpy as np
    import pandas as pd
    import scanpy as sc
    sc.settings.n_jobs = n_jobs  ## default = 8 
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    ####
    from sklearn.metrics import silhouette_score
    from sklearn.metrics import silhouette_samples
    #############################################################################   Import everything needed END
    
    
    ################################################################################# directory set up 
    os.makedirs(output_dir+output_prefix, exist_ok=True) ############################This is where all files related to the dataset associated with <output_prefix> will go
    
    
    os.makedirs(output_dir+output_prefix+'/tables/', exist_ok=True) ############################This is where table of sillhoutte scores will go (not used in this function)
    dataset_tables_output_directory=output_dir+output_prefix+'/tables/'  
    
    os.makedirs(output_dir+output_prefix+'/figures/', exist_ok=True)  ############################This is where  all figures related to the dataset associated with <output_prefix> will go
    dataset_figures_output_directory=output_dir+output_prefix+'/figures/'

    sc.settings.figdir=dataset_figures_output_directory  ############################This sets the scanpy figure output directory  
    
    ################################################################################# directory set up END   
    
    
    #### ################################################################## Leiden based clustering and sillhouette scoreing 
    ##################### Leiden based clustering
    # leiden resolution set in parameters
    print(f'leiden_res parameter set to {leiden_res}')   # Print the leiden resolution parameter used  
    sc.tl.leiden(adata,resolution=leiden_res)  # perofrom the leiden clustering 
    ##################### Leiden based clustering ##### END




    ##################### sillhouette scoreing
    samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden'])
    adata.obs['silhoutte']=samples_silhoutte_scores.tolist()
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    ##################### sillhouette scoreing  #####END
    ####################################################################### Leiden based clustering and sillhouette scoreing END
    
    ####################################################################### MAKE PLOTS
    
    ###################### umap and sillhouette scoreing graph results of final leiden resolution setting 
    fig_leiden_res_cluster_scores, (UMAP_final,ax_final,UMAP_sil,pca_leiden) = plt.subplots(nrows=1, ncols=4, figsize=(20,5), gridspec_kw={'wspace':0.4})

    ######################################################################################## plot  1) Umap plot colored by the leiden clusters.
    UMAP_final=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata), ax=UMAP_final,palette=sc.pl.palettes.vega_20_scanpy ,
                          show=False)

    cluster_silhouette_score_list=[]
    for i in range(0, cluster_number):
        cluster_silhouette_score=adata.obs['silhoutte'].loc[adata.obs['leiden']==str(i)].mean()
        cluster_silhouette_score_list.append(cluster_silhouette_score) 

    #  cluster scores
    pre_scores=cluster_silhouette_score_list
    pre_y_pos = np.arange(len(cluster_silhouette_score_list))
    ax_final.barh(pre_y_pos,pre_scores)
    ax_final.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_final.invert_yaxis()  # labels read top-to-bottom
    ax_final.set_title('Cluster Scores')  ############################################## plot   2) Bar graph of sillhouette score of each cluster 

    # #################################################################################plot 3) Umap plot colored by sillhouette score of each data point (cell) 
    UMAP_sil=sc.pl.umap(adata, color=['silhoutte'], ax=UMAP_sil,#palette=sc.pl.palettes.vega_20_scanpy,
                        show=False) 
    
    pca_leiden=sc.pl.pca(adata, color='leiden', ax=pca_leiden, show=False)  ########## plot 4) PCA plot of first two PCAs colored by leiden clusters 
    
    if savefig ==True: ############################ saves plot if true 
        os.makedirs(output_dir+output_prefix+'/', exist_ok=True)
        fig_leiden_res_cluster_scores.savefig(dataset_figures_output_directory+'leiden_res_'+str(leiden_res)+'_cluster_scores.pdf')
    
    ####################################################################### MAKE PLOTS END
    
    return adata # return adata objest with leiden clustering added to adata.obs['leiden']



""
def silhouette_walk_Largest_drop(
    adata: anndata.AnnData | None = None,
    leiden_res_start: float = 0.025,
    leiden_res_end: float = 2.0,
    leiden_res_step: float = 0.025,
    print_per_step: bool = False,
    n_jobs: int = 1,
    savetable: bool = False,
    savefig: bool = False,
    output_dir: str = "./figures/",
    output_prefix: str = "adata",
):
    """
    defualt values :
    leiden_res_start=0.025,leiden_res_end=2,leiden_res_step=0.025,print_per_step=False,
    n_jobs=8,savetable=False,savefig=False,output_dir="./figures/",output_prefix="adata"
    Thsi function walks through a range of leiden clustering resolutions defined by leiden_res_start,leiden_res_end,leiden_res_step.
    At step through walk 4 values are added to a table
    1. leiden each resoltion value used for clutering at step n
    2.the number of clusters at step n
    3.average sillhoutte score for all observations at step n
    4.average sillhoutte score of the observations in each cluster at step n
    After the walk is complete the  delta_silhouette_score is caluclated for each row of the table (skiping first row [1:])
    This difference between the average sillhoutte score for all observations at step n-1 minus the average sillhoutte score for all observations at step n
    the delta_silhouette_score is added as the 5th column of the table 
    The leiden resolution values before and after the largest drop in average sillhoutte score for all observations are indentified 
    For this pair of leiden resolutions leiden clustering is perofrmed and  The umaps are the plotted along with a bar chart of average sillhoutte score of the observations in each cluster
    """
        
    import os
    import numpy as np
    import pandas as pd
    import scanpy as sc

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    ####
    from sklearn.metrics import silhouette_score
    from sklearn.metrics import silhouette_samples
    ##########


    ####
    sc.settings.verbosity = 1 
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    sc.settings.n_jobs = n_jobs  ## 

    ################################################################################# directory set up 
    os.makedirs(output_dir+output_prefix, exist_ok=True)
    
    os.makedirs(output_dir+output_prefix+'/tables/', exist_ok=True)
    dataset_tables_output_directory=output_dir+output_prefix+'/tables/'
    
    os.makedirs(output_dir+output_prefix+'/figures/', exist_ok=True)
    dataset_figures_output_directory=output_dir+output_prefix+'/figures/'

    sc.settings.figdir=dataset_figures_output_directory
    ################################################################################# directory set up END   

    ###### cell checks the cluster at before and after  largest drops in average Sil score
    # leiden resolutioni set at value before the pargest drop  if not changed after


    silhouette_score_results=pd.DataFrame()
    leiden_res_range= np.arange(leiden_res_start,leiden_res_end, leiden_res_step).tolist()
    leiden_resolution_list=[]
    cluster_number_list=[]
    avg_silhouette_score_list=[]
    cluster_silhouette_score_list_all=[]
    for leiden_resolution in leiden_res_range:
        sc.tl.leiden(adata,resolution=leiden_resolution)
        cluster_number=len(set(adata.obs['leiden'].tolist()))
        samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden'])
        adata.obs['silhoutte']=samples_silhoutte_scores.tolist()
        average_silhouette_score=adata.obs['silhoutte'].mean()
        cluster_silhouette_score_list=[] 
        for i in range(0, cluster_number):
            cluster_silhouette_score=adata.obs['silhoutte'].loc[adata.obs['leiden']==str(i)].mean()
            cluster_silhouette_score_list.append(cluster_silhouette_score) 
        # add new values for this leiden resolution i into lists (to form datafrfame columns)
        leiden_resolution_list.append(leiden_resolution)
        cluster_number_list.append(cluster_number)
        avg_silhouette_score_list.append(average_silhouette_score)
        cluster_silhouette_score_list_all.append(cluster_silhouette_score_list)
        if print_per_step==True:         
            print(f' Average silhoutte score = {average_silhouette_score} for {cluster_number} clusters at leiden resolution of {leiden_res} {os.linesep} list of cluster scores = {cluster_silhouette_score_list}')

    # Make data frame from lists
    silhouette_score_results['leiden_resolution']=leiden_resolution_list
    silhouette_score_results['cluster_number']=cluster_number_list
    silhouette_score_results['silhouette_score']=avg_silhouette_score_list
    silhouette_score_results['cluster_silhouette_scores']=cluster_silhouette_score_list_all

    # make column with the difference in avg silhouette score between current and last leident resolution
    delta_silhouette_score_list=[0]
    for i in silhouette_score_results.index.tolist()[1:]:
        delta_silhouette_score=silhouette_score_results["silhouette_score"][i-1]-silhouette_score_results["silhouette_score"][i]
        delta_silhouette_score_list.append(delta_silhouette_score)
    silhouette_score_results['delta_silhouette_score']=delta_silhouette_score_list

    # amek data frame with sorted by delta_silhouette_score
    sil_results_maxD_sort=silhouette_score_results.sort_values('delta_silhouette_score',ascending=False)


    ##################### Leiden res values and cluster sil scores for pre and post LARGEST drop in avg sil score
    post_drop_1=sil_results_maxD_sort.index.tolist()[0] ### index values of LARGEST drop in sil score
    ### leiden resolutions beofore and after Biggest drop in sil score
    leiden_res_pre_drop1=float(silhouette_score_results.iloc[[post_drop_1-1]]['leiden_resolution'])
    leiden_res_post_drop1=float(silhouette_score_results.iloc[[post_drop_1]]['leiden_resolution'])
    #### cluster_silhouette_scores beofore and after Biggest drop in sil score
    cluster_silhouette_scores_pre1=pd.DataFrame(silhouette_score_results.iloc[[post_drop_1-1]]['cluster_silhouette_scores'].tolist()[-1])
    cluster_silhouette_scores_post1=pd.DataFrame(silhouette_score_results.iloc[[post_drop_1]]['cluster_silhouette_scores'].tolist()[-1])


    ######################################### graph results 
    # create objects
    fig = plt.figure(constrained_layout=True, figsize=(20,30),)
    gs = GridSpec(9, 4, figure=fig, )

    ###############
    ax= fig.add_subplot(gs[0, :])
    UMAP_pre = fig.add_subplot(gs[1, 0])
    ax_pre = fig.add_subplot(gs[1, 1])
    UMAP_post = fig.add_subplot(gs[1, 2])
    ax_post = fig.add_subplot(gs[1, 3])
    
    

    ################# graph average sil core vs leiden resolution 
    ax.plot(silhouette_score_results.leiden_resolution, silhouette_score_results.cluster_number, color="red", marker="o")
    ax.set_xlabel("leiden_resolution",fontsize=14)
    ax.set_ylabel("cluster_number",color="red",fontsize=14)
    ax2=ax.twinx()
    ax2.plot(silhouette_score_results.leiden_resolution, silhouette_score_results["silhouette_score"],color="blue",marker="o")
    ax2.set_ylabel("silhouette_score",color="blue",fontsize=14)

    ############################# LARGEST drop 

    # post drop #1 UMAP
    leiden_res=leiden_res_post_drop1
    print(f'leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_post=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),ax=UMAP_post,palette=sc.pl.palettes.vega_20_scanpy,
                         show=False)
    # post drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_post1
    post_scores=cluster_silhouette_scores_N[0].tolist()
    post_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_post.barh(post_y_pos, post_scores)
    ax_post.set_yticks(post_y_pos, labels=post_y_pos)
    ax_post.invert_yaxis()  # labels read top-to-bottom
    ax_post.set_title('Cluster Scores Post Drop1')

    # pre drop #1 UMAP ### this leaves leiden clustering set at value just beofre largest drop 
    leiden_res=leiden_res_pre_drop1
    print(f'Final leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_pre=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                        ax=UMAP_pre, show=False)
    # pre drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_pre1
    pre_scores=cluster_silhouette_scores_N[0].tolist()
    pre_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_pre.barh(pre_y_pos,pre_scores)
    ax_pre.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_pre.invert_yaxis()  # labels read top-to-bottom
    ax_pre.set_title('Cluster Scores Pre Drop1')

    #############################

    ######################################### grpah results END
    
    
    ############################# Save Table 
    if savetable ==True:
       # os.makedirs(output_dir+output_prefix+'/', exist_ok=True)
        silhouette_score_results.to_csv((dataset_tables_output_directory+output_prefix+'silhouette_score_results.csv'))
    ############################# Save Table END
    
    ############################# Save figure 
    if savefig ==True:
       # os.makedirs(output_dir+output_prefix+'/', exist_ok=True)
        fig.savefig(dataset_figures_output_directory+output_prefix+'_sil_walk_drop_1.pdf')
    ############################# Save figure END 
    
    
    
    return adata
    




""
def silhouette_walk_4_Largest_drops(
    adata: anndata.AnnData | None = None,
    leiden_res_start: float = 0.025,
    leiden_res_end: float = 2.0,
    leiden_res_step: float = 0.025,
    print_per_step: bool = False,
    n_jobs: int = 1,
    savetable: bool = False,
    savefig: bool = False,
    output_dir: str = "./figures/",
    output_prefix: str = "adata"
):
    """
    defualt values :
    leiden_res_start=0.025,leiden_res_end=2,leiden_res_step=0.025,print_per_step=False,
    n_jobs=8,savetable=False,savefig=False,output_dir="./figures/",output_prefix="adata"
    Thsi function walks through a range of leiden clustering resolutions defined by leiden_res_start,leiden_res_end,leiden_res_step.
    At step through walk 4 values are added to a table
    1. leiden each resoltion value used for clutering at step n
    2.the number of clusters at step n
    3.average sillhoutte score for all observations at step n
    4.average sillhoutte score of the observations in each cluster at step n
    After the walk is complete the  delta_silhouette_score is caluclated for each row of the table (skiping first row [1:])
    This difference between the average sillhoutte score for all observations at step n-1 minus the average sillhoutte score for all observations at step n
    the delta_silhouette_score is added as the 5th column of the table 
    The leiden resolution values before and after the the 4 largest drops in average sillhoutte score for all observations are indentified 
    For each of these 4 pairs of leiden resolutions  leiden clustering is perofrmed and  The umaps are the plotted along with a bar chart of average sillhoutte score of the observations in each cluster
    """
    import os
    import numpy as np
    import pandas as pd
    import scanpy as sc

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    ####
    from sklearn.metrics import silhouette_score
    from sklearn.metrics import silhouette_samples
    ##########


    ####
    sc.settings.verbosity = 1 
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')
    sc.settings.n_jobs = n_jobs  ## try this sometime

    ################################################################################# directory set up 
    if savefig:
        os.makedirs(output_dir+output_prefix, exist_ok=True)

        os.makedirs(output_dir+output_prefix+'/tables/', exist_ok=True)
        dataset_tables_output_directory=output_dir+output_prefix+'/tables/'

        os.makedirs(output_dir+output_prefix+'/figures/', exist_ok=True)
        dataset_figures_output_directory=output_dir+output_prefix+'/figures/'

        sc.settings.figdir=dataset_figures_output_directory
    ################################################################################# directory set up END   


###### cell checks the cluster at the 4 largest drops in average Silscore
# leiden resolutioni set at value before the pargest drop  if not changed after


    leiden_res_range= np.arange(leiden_res_start,leiden_res_end, leiden_res_step).tolist()

    silhouette_score_results=pd.DataFrame()
    leiden_resolution_list=[]
    cluster_number_list=[]
    avg_silhouette_score_list=[]
    cluster_silhouette_score_list_all=[]

    for leiden_resolution in leiden_res_range:
        sc.tl.leiden(adata,resolution=leiden_resolution)
        cluster_number=len(set(adata.obs['leiden'].tolist()))
        samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden'])
        adata.obs['silhoutte']=samples_silhoutte_scores.tolist()
        average_silhouette_score=adata.obs['silhoutte'].mean()
        cluster_silhouette_score_list=[] 
        for i in range(0, cluster_number):
            cluster_silhouette_score=adata.obs['silhoutte'].loc[adata.obs['leiden']==str(i)].mean()
            cluster_silhouette_score_list.append(cluster_silhouette_score) 
        #add new values for this leiden resolution i into lists (to form datafrfame columns)
        leiden_resolution_list.append(leiden_resolution)
        cluster_number_list.append(cluster_number)
        avg_silhouette_score_list.append(average_silhouette_score)
        cluster_silhouette_score_list_all.append(cluster_silhouette_score_list)
        if print_per_step==True:         
            print(f' Average silhoutte score = {average_silhouette_score} for {cluster_number} clusters at leiden resolution of {leiden_res} {os.linesep} list of cluster scores = {cluster_silhouette_score_list}')

    # Make data frame from lists
    silhouette_score_results['leiden_resolution']=leiden_resolution_list
    silhouette_score_results['cluster_number']=cluster_number_list
    silhouette_score_results['silhouette_score']=avg_silhouette_score_list
    silhouette_score_results['cluster_silhouette_scores']=cluster_silhouette_score_list_all

    # make column with the difference in avg silhouette score between current and last leident resolution
    delta_silhouette_score_list=[0]
    for i in silhouette_score_results.index.tolist()[1:]:
        delta_silhouette_score=silhouette_score_results["silhouette_score"][i-1]-silhouette_score_results["silhouette_score"][i]
        delta_silhouette_score_list.append(delta_silhouette_score)
    silhouette_score_results['delta_silhouette_score']=delta_silhouette_score_list

    # amek data frame with sorted by delta_silhouette_score
    sil_results_maxD_sort=silhouette_score_results.sort_values('delta_silhouette_score',ascending=False)


    ##################### Leiden res values and cluster sil scores for pre and post LARGEST drop in avg sil score
    post_drop_1=sil_results_maxD_sort.index.tolist()[0] ### index values of LARGEST drop in sil score
    ### leiden resolutions beofore and after Biggest drop in sil score
    leiden_res_pre_drop1=float(silhouette_score_results.iloc[[post_drop_1-1]]['leiden_resolution'])
    leiden_res_post_drop1=float(silhouette_score_results.iloc[[post_drop_1]]['leiden_resolution'])
    #### cluster_silhouette_scores beofore and after Biggest drop in sil score
    cluster_silhouette_scores_pre1=pd.DataFrame(silhouette_score_results.iloc[[post_drop_1-1]]['cluster_silhouette_scores'].tolist()[-1])
    cluster_silhouette_scores_post1=pd.DataFrame(silhouette_score_results.iloc[[post_drop_1]]['cluster_silhouette_scores'].tolist()[-1])

    ##################### Leiden res values and cluster sil scores for pre and post SECOND largest drop in avg sil score
    post_drop_2=sil_results_maxD_sort.index.tolist()[1] ### index values of LARGEST drop in sil score
    ### leiden resolutions beofore and after Biggest drop in sil score
    leiden_res_pre_drop2=float(silhouette_score_results.iloc[[post_drop_2-1]]['leiden_resolution'])
    leiden_res_post_drop2=float(silhouette_score_results.iloc[[post_drop_2]]['leiden_resolution'])
    #### cluster_silhouette_scores beofore and after Biggest drop in sil score
    cluster_silhouette_scores_pre2=pd.DataFrame(silhouette_score_results.iloc[[post_drop_2-1]]['cluster_silhouette_scores'].tolist()[-1])
    cluster_silhouette_scores_post2=pd.DataFrame(silhouette_score_results.iloc[[post_drop_2]]['cluster_silhouette_scores'].tolist()[-1])

    ##################### Leiden res values and cluster sil scores for pre and post THIRD largest drop in avg sil score
    post_drop_3=sil_results_maxD_sort.index.tolist()[2] ### index values of LARGEST drop in sil score
    ### leiden resolutions beofore and after Biggest drop in sil score
    leiden_res_pre_drop3=float(silhouette_score_results.iloc[[post_drop_3-1]]['leiden_resolution'])
    leiden_res_post_drop3=float(silhouette_score_results.iloc[[post_drop_3]]['leiden_resolution'])
    #### cluster_silhouette_scores beofore and after Biggest drop in sil score
    cluster_silhouette_scores_pre3=pd.DataFrame(silhouette_score_results.iloc[[post_drop_3-1]]['cluster_silhouette_scores'].tolist()[-1])
    cluster_silhouette_scores_post3=pd.DataFrame(silhouette_score_results.iloc[[post_drop_3]]['cluster_silhouette_scores'].tolist()[-1])

    ##################### Leiden res values and cluster sil scores for pre and FOURTH largest drop in avg sil score
    post_drop_4=sil_results_maxD_sort.index.tolist()[3] ### index values of LARGEST drop in sil score
    ###leiden resolutions beofore and after Biggest drop in sil score
    leiden_res_pre_drop4=float(silhouette_score_results.iloc[[post_drop_4-1]]['leiden_resolution'])
    leiden_res_post_drop4=float(silhouette_score_results.iloc[[post_drop_4]]['leiden_resolution'])
    #### cluster_silhouette_scores beofore and after Biggest drop in sil score
    cluster_silhouette_scores_pre4=pd.DataFrame(silhouette_score_results.iloc[[post_drop_4-1]]['cluster_silhouette_scores'].tolist()[-1])
    cluster_silhouette_scores_post4=pd.DataFrame(silhouette_score_results.iloc[[post_drop_4]]['cluster_silhouette_scores'].tolist()[-1])



    ######################################### grpah results 
    # create objects
    fig = plt.figure(constrained_layout=True, figsize=(20,30),)
    
    
    gs = GridSpec(9, 4, figure=fig, )

    ax= fig.add_subplot(gs[0, :])
    UMAP_pre = fig.add_subplot(gs[1, 0])
    ax_pre = fig.add_subplot(gs[1, 1])
    UMAP_post = fig.add_subplot(gs[1, 2])
    ax_post = fig.add_subplot(gs[1, 3])
  


    #############################SECOND largest drop 
    UMAP_pre2 = fig.add_subplot(gs[2, 0])
    ax_pre2 = fig.add_subplot(gs[2, 1])
    UMAP_post2 = fig.add_subplot(gs[2, 2])
    ax_post2 = fig.add_subplot(gs[2, 3])
    #############################THIRD largest drop 
    UMAP_pre3 = fig.add_subplot(gs[3, 0])
    ax_pre3 = fig.add_subplot(gs[3, 1])
    UMAP_post3 = fig.add_subplot(gs[3, 2])
    ax_post3 = fig.add_subplot(gs[3, 3])
    #############################FOURTH largest drop
    UMAP_pre4 = fig.add_subplot(gs[4, 0])
    ax_pre4 = fig.add_subplot(gs[4, 1])
    UMAP_post4 = fig.add_subplot(gs[4, 2])
    ax_post4 = fig.add_subplot(gs[4, 3])


    ax.plot(silhouette_score_results.leiden_resolution, silhouette_score_results.cluster_number, color="red", marker="o")
    ax.set_xlabel("leiden_resolution",fontsize=14)
    ax.set_ylabel("cluster_number",color="red",fontsize=14)
    ax2=ax.twinx()
    ax2.plot(silhouette_score_results.leiden_resolution, silhouette_score_results["silhouette_score"],color="blue",marker="o")
    ax2.set_ylabel("silhouette_score",color="blue",fontsize=14)

    ############################# LARGEST drop 
    # pre drop #1 UMAP
    leiden_res=leiden_res_pre_drop1
    print(f' PRE drop 1 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_resolution}')
    UMAP_pre=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                        ax=UMAP_pre, show=False)
    # pre drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_pre1
    pre_scores=cluster_silhouette_scores_N[0].tolist()
    pre_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_pre.barh(pre_y_pos,pre_scores)
    ax_pre.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_pre.invert_yaxis()  # labels read top-to-bottom
    ax_pre.set_title('Cluster Scores Pre Drop1')
    # post drop #1 UMAP
    leiden_res=leiden_res_post_drop1
    print(f'POST drop 1 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_resolution}')
    UMAP_post=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                         ax=UMAP_post, show=False)
    # post drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_post1
    post_scores=cluster_silhouette_scores_N[0].tolist()
    post_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_post.barh(post_y_pos, post_scores)
    ax_post.set_yticks(post_y_pos, labels=post_y_pos)
    ax_post.invert_yaxis()  # labels read top-to-bottom
    ax_post.set_title('Cluster Scores Post Drop1')
    #############################



    ############################# SECOND largest drop 
    # pre drop #1 UMAP
    leiden_res=leiden_res_pre_drop2
    print(f' PRE drop 2 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_pre2=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy, 
                         ax=UMAP_pre2, show=False)
    # pre drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_pre2
    pre_scores=cluster_silhouette_scores_N[0].tolist()
    pre_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_pre2.barh(pre_y_pos,pre_scores)
    ax_pre2.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_pre2.invert_yaxis()  # labels read top-to-bottom
    ax_pre2.set_title('Cluster Scores Pre Drop2')
    # post drop #1 UMAP
    leiden_res=leiden_res_post_drop2
    print(f'POST drop 2 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_post2=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                          ax=UMAP_post2, show=False)
    # post drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_post2
    post_scores=cluster_silhouette_scores_N[0].tolist()
    post_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_post2.barh(post_y_pos, post_scores)
    ax_post2.set_yticks(post_y_pos, labels=post_y_pos)
    ax_post2.invert_yaxis()  # labels read top-to-bottom
    ax_post2.set_title('Cluster Scores Post Drop2')
    #############################


    ############################# THIRD largest drop 
    # pre drop #1 UMAP
    leiden_res=leiden_res_pre_drop3
    print(f' PRE drop 3 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_pre3=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                         ax=UMAP_pre3, show=False)
    # pre drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_pre3
    pre_scores=cluster_silhouette_scores_N[0].tolist()
    pre_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_pre3.barh(pre_y_pos,pre_scores)
    ax_pre3.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_pre3.invert_yaxis()  # labels read top-to-bottom
    ax_pre3.set_title('Cluster Scores Pre Drop3')
    # post drop #1 UMAP
    leiden_res=leiden_res_post_drop3
    print(f'POST drop 3 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_post3=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                          ax=UMAP_post3, show=False)
    # post drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_post3
    post_scores=cluster_silhouette_scores_N[0].tolist()
    post_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_post3.barh(post_y_pos, post_scores)
    ax_post3.set_yticks(post_y_pos, labels=post_y_pos)
    ax_post3.invert_yaxis()  # labels read top-to-bottom
    ax_post3.set_title('Cluster Scores Post Drop3')
    #############################


    ############################# FOURTH largest drop
    # pre drop #1 UMAP
    leiden_res=leiden_res_pre_drop4
    print(f'PRE drop 4 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_pre4=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy, 
                         ax=UMAP_pre4, show=False)
    # pre drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_pre4
    pre_scores=cluster_silhouette_scores_N[0].tolist()
    pre_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_pre4.barh(pre_y_pos,pre_scores)
    ax_pre4.set_yticks(pre_y_pos, labels=pre_y_pos)
    ax_pre4.invert_yaxis()  # labels read top-to-bottom
    ax_pre4.set_title('Cluster Scores Pre Drop4')
    # post drop #1 UMAP
    leiden_res=leiden_res_post_drop4
    print(f' POST drop 4 leiden_res parameter set to {leiden_res}')
    sc.tl.leiden(adata,resolution=leiden_res)
    cluster_number=len(set(adata.obs['leiden'].tolist()))
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f' Average silhoutte score = {silhouette_score_adata} for {cluster_number} clusters at leiden resolution of {leiden_res}')
    UMAP_post4=sc.pl.umap(adata, color='leiden',title='Leiden.res='+str(leiden_res)+' Avg.sil.='+str(silhouette_score_adata),palette=sc.pl.palettes.vega_20_scanpy,
                          ax=UMAP_post4, show=False)
    # post drop #1 cluster scores
    cluster_silhouette_scores_N=cluster_silhouette_scores_post4
    post_scores=cluster_silhouette_scores_N[0].tolist()
    post_y_pos = np.arange(len(cluster_silhouette_scores_N[0].tolist()))
    ax_post4.barh(post_y_pos, post_scores)
    ax_post4.set_yticks(post_y_pos, labels=post_y_pos)
    ax_post4.invert_yaxis()  # labels read top-to-bottom
    ax_post4.set_title('Cluster Scores Post Drop4')
    #############################
    ######################################### graph results END


    print(f'Current leiden_res parameter set to value before LARGEST (drop1) {leiden_res_pre_drop1}')
    sc.tl.leiden(adata,resolution=leiden_res_pre_drop1)
    silhouette_score_adata=silhouette_score(adata.obsm['X_pca'], adata.obs['leiden'],)
    print(f'Current avg. silhouette_score is  {silhouette_score_adata}')
    
    
    ############################# Save Table 
    if savetable ==True:
       # os.makedirs(output_dir+output_prefix+'/', exist_ok=True)
        silhouette_score_results.to_csv((dataset_tables_output_directory+output_prefix+'silhouette_score_results.csv'))
    ############################# Save Table END
    
    ############################# Save Figure 
    if savefig ==True:
        #os.makedirs(output_dir+output_prefix+'/', exist_ok=True)
        fig.savefig(dataset_figures_output_directory+output_prefix+'_sil_walk_drop_4.pdf')
    ############################# Save Figure end
    
    return adata



# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]