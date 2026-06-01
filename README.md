# single_cell_python_tools

## Documentation

- Published docs: [gitbenlewis.github.io/single_cell_python_tools](https://gitbenlewis.github.io/single_cell_python_tools/)
- Source docs: [docs/README.md](docs/README.md)
- Start with [installation](docs/installation.md), [quickstart](docs/quickstart.md), [DATASET_class](docs/DATASET_class.md), and the [API reference](docs/api-reference.md).

Mostly scanpy wrappers and tools that use adata objects to stream line single cell analysis 

# Installation
      git clone  <this repo>
      cd <this repo>
      conda env remove -n sctl  # Removes any existing 'sctl' environment
      conda env create -f sctl.yaml
      conda activate sctl

      or 
      
      git clone  <this repo>
      cd <this repo>
      pip install -e . # Installs the package in editable mode for development

# Environment update
      conda env update -f sctl.yaml

# After installation, verify the package:
    python -c "import single_cell_python_tools as sctl; print(sctl.__version__)"

# Package import
      import single_cell_python_tools as sctl
      print(sctl.DATASET_class.__doc__) # access DATASET_class via sctl.DATASET_class
      from single_cell_python_tools.DATASET_class import *
      print(DATASET_class.__doc__) # after to current namespace access DATASET_class via DATASET_class
   
# Example_notebooks/
 (most usefull starting point)
check out the example notebook folders
## /01_PBMC3k: QC,preprocess,clustering, visuallzation, DEG , GSEA
###  00_sctl_functions_preprocessing_PBMC3K.ipynb 
   #### functions operate a separate adata object
      - sctl.function(adata,**adata.uns["parameters"]) 
      Ex:
      #  data download handled elsewhere
      adata = sc.read_10x_mtx(data_directory_path)
      adata.uns["parameters"]=pbmc3k_parameters
      sctl.pp.basic_filitering(adata,**adata.uns["parameters"])
      sctl.pp.annotate_n_view_adata_raw_counts(adata,**adata.uns["parameters"])
      sctl.pp.filter_cells_by_anotated_QC_gene(adata,**adata.uns["parameters"])
      sctl.pp.remove_genes(adata,**adata.uns["parameters"])
      sctl.pp.process2scaledPCA(adata,**adata.uns["parameters"])
      sctl.pp.leiden_clustering(adata,**adata.uns["parameters"])
      sctl.pp.silhouette_score_n_plot(adata,**adata.uns["parameters"])
      import scanpy as sc
      sc.pl.umap(adata, color=adata.uns["parameters"]['umap_marker_gene_list'])

### 01_sctl_DATASET_class_preprocessing_PBMC3K.ipynb
#### class object for each dataset, methods for analysis 
      - sctl_DATASET_class().basic_filitering().clustering().etc()
      Full pipeline Ex:
      pbmc3k_sctl_DATASET=sctl.DATASET_class(parameters=pbmc3k_parameters)
      (
      pbmc3k_sctl_DATASET
      # datadown load
      .download_data()
      .unpack_tar()
      .load_data()
      # QC data
      .basic_filitering()
      .annotate_n_view_adata_raw_counts()
      .filter_cells_by_anotated_QC_gene()
      .remove_genes()
      .annotate_n_view_adata_raw_counts()
      ### transform and normalize
      .process2scaledPCA()
      ### clustering
      .leiden_clustering()
      .rename_leiden_clusters()
      .silhouette_score_n_plot()
      ### UMAP plots
      .marker_gene_umaps()
      )
### Example parameters dictionary
      # can store dictionary in adata.uns
      bmc3k_parameters={
      # meta
      "download_url":'http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz',
      "download_output_dir":'./data/',
      "download_output_filename": None,
      "input_file_path":'./data/filtered_gene_bc_matrices/hg19',
      'file_format':'10x',
      'dataset_prefix_for_10x_triplets':'',# include all characters up to barcodes/features (underscore)
      "output_prefix":'EX_pbmc3k_raw_sctl_DATASET_',
      "output_dir":'./write/',
      'make_empty_output_dirs': False,  # make empty output dirs if True or if output_dir is not none
      "n_jobs": 1,
      'organism' : 'human',

      #### dataset specfic parameters
      ###Basic filters
      "filter_genes_min_cells":3,  # min of cells a gene is detected in else gene is tossed out default 3
      "filter_genes_min_counts":0, # min  of counts a gene must have to pass basic filter default 0
      "filter_cells_min_genes":200, # min  of genes detected or else cell/observation is tossed out default 200
      "filter_cells_min_counts":0, # min  of counts detected or else cell/observation is tossed out default 0  
         
      ####Filter  on off switches
      "filter_ncount" : True,
      "filter_pct_mt" : True,
      "filter_pct_ribo" : False,
      "filter_pct_hb" : False,
      "filter_pct_malat1":False,
      "filter_HVG" : True,

      ###less than filter percent 
      "n_genes_bycounts" : 2500, #less than filter
      "percent_mt" : 5, #less than filter
      "percent_ribo" : 100, #less than filter
      "percent_malat1": 100, #less than filter
      "percent_hb" : 100,  #less than filter

      ###Greater than filter percent
      "over_n_genes_bycounts" : 200, #greater than filter
      "over_percent_mt" : 0, #greater than filter
      "over_percent_ribo" : 0, #greater than filter
      "over_percent_malat1": 0, #greater than filter
      "over_percent_hb" : 0 , #greater than filter

      ###Remove gene sets  on off switches
      "remove_MALAT1" : False, 
      "remove_MT" : False ,
      "remove_HB" : False,
      "remove_RP_SL" : False ,
      "remove_MRP_SL" : False,

      #### processing parameters and options
      "filter_genes_min_counts_normed":0,
      "normalize_total_target_sum" : 1e4,  # scanpy  default 1e4
      "HVG_min_mean"  :  0.0125, # scanpy  default 0.0125
      "HVG_max_mean"  :  3, # scanpy  default 3
      "HVG_min_disp"  :  0.5, # scanpy  default 0.5
      "logarithmize":True, # scanpy default True
      "scale":True, # scanpy default True
      "scale_max_std_value":10, # 10 often used
      "save_counts_layer": True,

      ####regression on off switches
      "regress_mt" : True,
      "regress_ribo" : False,
      "regress_malat1":False,
      "regress_hb" : False,
      "regress_cell_cycle_score" : False,

      ###clustering parameters for clusters
      "number_of_PC" : 40, ### dataset demensionality 
      "number_of_neighbors" : 10,
      "leiden_res" : .9, #leiden clustering resolution


      # UMAP graph parameters
      'umap_marker_gene':True,
      'umap_marker_gene_list': ['IL7R','CD14','LYZ', 'MS4A1','CD8A','GNLY','NKG7','FCGR3A','MS4A7','FCER1A','CST3','PPBP'],
      'ncols': 4,
      'vmax':None,
      #cluster naming parameters
      'rename_cluster': True,
      ####  'NK' and 'FCGR3A Monocytes' seem to switch places sometimes .. not sure if its random or a bug, NK cluster should be close to CD8 T cells
      #'new_cluster_names' : [  'CD4 T', 'CD14 Monocytes', 'B','CD8 T',   'FCGR3A Monocytes','NK','Dendritic', 'Megakaryocytes'],
      'new_cluster_names' : [  'CD4 T', 'CD14 Monocytes', 'B','CD8 T',   'NK','FCGR3A Monocytes','Dendritic', 'Megakaryocytes'],  
      }


# Getting started


## Code to add to the start of a jupyter notebook in order to run the package
``` 
# option 1 use an absolute path to direcotry containing python packages from github
repo_parent_dir='home/ubuntu/projects/github_repos'
# option 2 use a relative path that will work as long as a the direcotry containing python packages is two levels above the directory containing the jupyter notebook
repo_parent_dir='../../'

import sys
if repo_parent_dir not in sys.path:
    sys.path.append(repo_parent_dir)

# import the package
import single_cell_python_tools as sctl

``` 
## Directory structure

I recommend placing all your github repos (packages and analysis_project_repos) in the same parent directory that way you can use relative paths to access packages from your analysis_project_repos.\
This way the same notebook will work on your aws enviroment as as well as your local envirmoment ( without changing path variables).\
For example adding, adding '../../' to your sys.path from any of the notebooks in the example direcotry tree below will give you access to both the single_cell_python_tools and the second_favorite_package_repo_directory regardless of where the "github_repos" is located.

Recommended directory structure:
github_repos/
  ├── single_cell_python_tools/
  ├── second_favorite_package_repo/
  └── singlecell_analysis_project_repo/
      ├── singlecell_analysis_project_A/
      └── singlecell_analysis_project_B/



# single_cell_python_tools package functions

## Getting start will the package functions
### Check the doc strings of the functions for more info on each

`print(functionOrClass.__doc__) #(double underscores)`

`print(sctl.pp.basic_filitering.__doc__)`

`print(sctl.pp.annotate_n_view_adata_raw_counts.__doc__)`

`print(sctl.pp.filter_cells_by_anotated_QC_gene.__doc__)`

`print(sctl.pp.remove_genes.__doc__)`

`print(sctl.pp.process2scaledPCA.__doc__)`

`print(sctl.pp.leiden_clustering.__doc__)`

`print(sctl.pl.silhouette_score_n_plot.__doc__)`


### latest update is a modified version of scanpy's ingest funciton
`print(sctl.tl.ingest_verbose.__doc__)`

## Subpackage - single_cell_python_tools - preprocessing - sctl.pp.

`sctl.pp.basic_filitering()`
```
basic_filitering(adata,
                     filter_cells_min_counts=0,
                      filter_cells_min_genes=200,
                     filter_genes_min_cells=3,
                     filter_genes_min_counts=0,
                     **parameters):
```


`sctl.pp.annotate_QC_genes(adata)`

```
annotate the group of QC genes
### double HB gene annoation works.... maybe just give it a list of genes
#### code 
adata.var['mt'] = adata.var_names.str.startswith("MT-")  # mitochondrial genes as 'mt'
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL")) # ribosomal genes genes as 'ribo'
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)(S)]")) & ~adata.var_names.str.contains(("HBEGF")) 
# "^HB[^(P)" changed to "^HB[^(P)(S)" and  & ~adata_test.var_names.str.contains(("HBEGF")) added to remove HBS1L and HBEGF which are NOT memoglobin genes
adata.var['malat1'] = adata.var_names.str.contains(("MALAT1"))  # MALAT1 genes as 'malat1'
return adata
```


`sctl.pp.calculate_qc_metrics(adata)`
```
calculate_qc_metrics
# add code to check if genes already annotated 
#### code 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # mitocohndrial  genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True) # ribosomal genes
sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True) # hemoglobin genes.
sc.pp.calculate_qc_metrics(adata, qc_vars=['malat1'], percent_top=None, log1p=False, inplace=True) # MALAT1 gene.
return adata
```

`sctl.pp.plot_QC_metrics_scatter(adata)`
```
#### code 
figQC, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1 ,5,figsize=(20,4), gridspec_kw={'wspace':0.9})
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',ax=ax1, show=False) # plot number of dected genes vs total counts 
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',ax=ax2, show=False) #percent mt counts vs total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_ribo',ax=ax3, show=False) #percent ribo counts vs total counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_malat1',ax=ax4, show=False) #percent HB counts vs total count
sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb',ax=ax5, show=False) #percent HB counts vs total counts
```

`sctl.pp.plot_QC_metrics_violin(adata)`
```
#### code 
fig1, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(1 ,6,figsize=(20,4), gridspec_kw={'wspace':0.9})
sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4,ax=ax1, show=False)
sc.pl.violin(adata, ['total_counts'], jitter=0.4 ,ax=ax2, show=False)
sc.pl.violin(adata, [ 'pct_counts_mt'], jitter=0.4,ax=ax3, show=False) # mitocohndrial  genes
sc.pl.violin(adata, [ 'pct_counts_ribo'], jitter=0.4,ax=ax4, show=False) # ribosomal genes
sc.pl.violin(adata, [ 'pct_counts_malat1'], jitter=0.4,ax=ax5, show=False) # hemoglobin genes.
sc.pl.violin(adata, [ 'pct_counts_hb'], jitter=0.4,ax=ax6, show=False) # hemoglobin genes.  
```

`sctl.pp.plot_qc_metrics(adata)`
```
plot_qc_metrics of Annotated technical gene groups  and top 20 highly expressed
#### code 
plot_QC_metrics_violin(adata)  
plot_QC_metrics_scatter(adata) 
sc.pl.highest_expr_genes(adata, n_top=20, )
```

`sctl.pp.annotate_n_view_adata_raw_counts(adata)`
```
Annotate technical gene groups  and calculate qc metrics
#### code 
annotate_QC_genes(adata)
calculate_qc_metrics(adata)
plot_qc_metrics(adata)
```


`sctl.pp.filter_cells_by_anotated_QC_gene()`
```
filter_cells_by_anotated_QC_gene(adata,
                                    filter_ncount=True,
                                    n_genes_bycounts=10000,
                                    filter_pct_mt=True,
                                    percent_mt=20,
                                    over_percent_mt=0,
                                    filter_pct_ribo=False,
                                    percent_ribo=100,
                                    over_percent_ribo=0,
                                    filter_pct_hb=False,
                                    percent_hb=100,
                                    over_percent_hb=0,
                                    filter_pct_malat1=False,
                                    percent_malat1=100,
                                    over_percent_malat1=0,
                                      **parameters
                                    )
```

`sctl.pp.remove_genes()`
```
remove_genes(adata,
                remove_MALAT1=False,
                remove_MT=False,
                remove_HB=False,
                remove_RP_SL=False,
                remove_MRP_SL=False,
                 **parameters
                )
```

`sctl.pp.process2scaledPCA()`
```
 process2scaledPCA(adata,
                   normalize_total_target_sum=1e4,logarithmize=True,
                   filter_HVG=False,HVG_min_mean=0.0125, HVG_max_mean=3, HVG_min_disp=3,
                   regress_mt=False,regress_ribo=False,regress_malat1=False,regress_hb=False,
                   scale=True,PCA=True,
                      cell_cycle_score=True,
                   regress_cell_cycle_score=False,HVG_flavor='seurat',HVG_n_top_genes=1500,
                      **parameters)
```

`sctl.pp.leiden_clustering()`
```
#### code
leiden_clustering(adata, number_of_neighbors=10,number_of_PC=40, leiden_res=1,key_added='leiden', **parameters)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=leiden_res,key_added=key_added)
sc.pl.umap(adata, color=[key_added] )
```




## Subpackage - single_cell_python_tools - tools - sctl
* analysis helper functions (not QC or preprocessing)
 

## Subpackage - single_cell_python_tools - plots - sctl
`sctl.pl.silhouette_score_n_plot(adata,leiden_res='unk',**parameters)`
```
sctl.pl.silhouette_score_n_plot(adata,parameters,leiden_res='unk'):
> assumes ledien clusteirng to subset cells
> uses X_pca for silhoutte_scores
samples_silhoutte_scores=silhouette_samples(adata.obsm['X_pca'], adata.obs['leiden']

```

`sctl.pl.plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count")`
```
sctl.pl.plot_batch_obs_key_of_obs_key2(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4),flavor="pct_count")
makes two side by side bar charts, each bar is a batch_obs_key='batch category and each bar is stacked and colored by obs_key2="leiden"
left bar chart is fraction on y -axis 
right  bar chart is obs/cell count on y -axis 
flavor="pct_count"  >>> both charts
flavor="pct"  >>> only pct chart
flavor="count"  >>> only count chart

```

`sctl.pl.plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4))`
```
sctl.pl.plot_percent_obs_key2_per_batch_obs_key(adata,savefig=False,output_dir='./project/',output_prefix="dataset_",batch_obs_key='batch',obs_key2="leiden",figsize=(10,4)
This produce one column of individual bar charts (one chart for each catagory in batch_obs_key='batch'") # batch_obs_key="sample_ID is good to use
each bar chart show percentage of cells in "batch" assigned to obs_key2="leiden"
```




# Stand alone module - sctl_gex class- single_cell_python_tools - sctl_gex - sctl.sctl_gex

`sctl_gex(adata,**parameters)`
* see the notebook 
   * sctl_gex_class_preprocessing_PBMC3k.ipynb # single class object for preprocessing


# more coming soon
