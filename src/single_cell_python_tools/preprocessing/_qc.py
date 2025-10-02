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


### annotate_QC_genes and calculate_qc_metrics functions
def annotate_QC_genes(
    adata: anndata.AnnData | None = None,
    organism: str = 'human',
    **parameters: Any
):
    """
    annotate_QC_genes
    """
    logger.info(f"Organism is set to {organism}")
    if organism == 'human' or None:
        logger.info(f"annotating human genes")
        #print ('Organism is human, annotating human genes')
        adata.var['mt'] = adata.var_names.str.startswith("MT-")  # mitochondrial genes as 'mt'
        adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL")) # ribosomal genes genes as 'ribo'
        adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)(S)]")) & ~adata.var_names.str.contains(("HBEGF")) 
        # "^HB[^(P)" changed to "^HB[^(P)(S)" and  & ~adata_test.var_names.str.contains(("HBEGF")) added to remove HBS1L and HBEGF which are NOT memoglobin genes
        adata.var['malat1'] = adata.var_names.str.contains(("MALAT1"))  # MALAT1 genes as 'malat1'
    if organism == 'mouse':
        logger.info(f" annotating mouse QC genes")
        #print ('Organism is mouse, annotating mouse genes')
        adata.var['mt'] = adata.var_names.str.startswith("mt-")  # mitochondrial genes as 'mt'
        adata.var['ribo'] = adata.var_names.str.startswith(("Rps","Rpl")) # ribosomal genes genes as 'ribo'
        adata.var['hb'] = adata.var_names.str.contains(("^Hb[^(P)(S)]")) & ~adata.var_names.str.contains(("Hbegf"))
        adata.var['malat1'] = adata.var_names.str.contains(("Malat1"))  # MALAT1 genes as 'malat1'
    return 

def calculate_qc_metrics(
    adata: anndata.AnnData | None = None,
    **parameters: Any
):
    """
    calculate_qc_metrics
    # add code to check if genes already annotated  
    """
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # mitocohndrial  genes
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['ribo'], percent_top=None, log1p=False, inplace=True) # ribosomal genes
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['hb'], percent_top=None, log1p=False, inplace=True) # hemoglobin genes.
    #sc.pp.calculate_qc_metrics(adata, qc_vars=['malat1'], percent_top=None, log1p=False, inplace=True) # MALAT1 gene.
    sc.pp.calculate_qc_metrics( adata, qc_vars=["mt", "ribo", "hb",'malat1'],percent_top=None, inplace=True, log1p=False) # mitocohndrial,ribosomal,hemoglobin,MALAT1  genes
    return 

### annotate_QC_genes and MD_calculate_qc_metrics functions END


### plot functions
def plot_QC_metrics_scatter(adata: anndata.AnnData | None = None,):
    '''
    '''
    figQC, (ax1,ax2,ax3,ax4,ax5) = plt.subplots(1 ,5,figsize=(20,4), gridspec_kw={'wspace':0.9})
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',ax=ax1, show=False) # plot number of dected genes vs total counts 
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',ax=ax2, show=False) #percent mt counts vs total counts
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_ribo',ax=ax3, show=False) #percent ribo counts vs total counts
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_malat1',ax=ax4, show=False) #percent HB counts vs total count
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb',ax=ax5, show=False) #percent HB counts vs total counts 
    return

def plot_QC_metrics_violin(adata: anndata.AnnData | None = None,):
    '''
    #### code 
    '''
    fig1, (ax1,ax2,ax3,ax4,ax5,ax6) = plt.subplots(1 ,6,figsize=(20,4), gridspec_kw={'wspace':0.9})
    sc.pl.violin(adata, ['n_genes_by_counts'], jitter=0.4,ax=ax1, show=False)
    sc.pl.violin(adata, ['total_counts'], jitter=0.4 ,ax=ax2, show=False)
    sc.pl.violin(adata, [ 'pct_counts_mt'], jitter=0.4,ax=ax3, show=False) # mitocohndrial  genes
    sc.pl.violin(adata, [ 'pct_counts_ribo'], jitter=0.4,ax=ax4, show=False) # ribosomal genes
    sc.pl.violin(adata, [ 'pct_counts_malat1'], jitter=0.4,ax=ax5, show=False) # hemoglobin genes.
    sc.pl.violin(adata, [ 'pct_counts_hb'], jitter=0.4,ax=ax6, show=False) # hemoglobin genes.
    return

def plot_qc_metrics(adata: anndata.AnnData | None = None, **parameters: Any):
    """
    plot_qc_metrics of Annotated technical gene groups  and top 20 highly expressed
    """
    plot_QC_metrics_violin(adata)  
    plot_QC_metrics_scatter(adata) 
    sc.pl.highest_expr_genes(adata, n_top=20, )
    return
    

### plot functions end 


### multi funcitons

def annotate_n_view_adata_raw_counts(
    adata: anndata.AnnData | None = None,
    **parameters: Any):
    """
    Annotate technical gene groups  and calculate qc metrics
    """
    annotate_QC_genes(adata,**parameters)
    calculate_qc_metrics(adata)
    plot_qc_metrics(adata) 
    return 

### multi funcitons END 



### filter funcitons 

def basic_filitering(
    adata: anndata.AnnData | None = None,
    filter_cells_min_counts: int = 0,
    filter_cells_min_genes: int = 200,
    filter_genes_min_cells: int = 3,
    filter_genes_min_counts: int = 0,
    **parameters: Any):
    """ Basic Filtering
  
    """
    logger.info(f" {adata.n_obs} observations BEFORE Basic Filtering (min_genes and min_counts)")
    #print(f" {adata.n_obs} observations BEFORE Basic Filtering")
    sc.pp.filter_cells(adata, min_genes=filter_cells_min_genes)  #min_genes=over_n_genes_bycounts
    logger.info(f" Filtering cells  sc.pp.filter_cells(adata, min_genes={filter_cells_min_genes})")
    logger.info(f' {adata.n_obs} observations AFTER  min_genes Filtering, min_genes={filter_cells_min_genes}')
    #print(f'Filtering cells pp.filter_cells(adata, min_genes=filter_cells_min_genes)  Cells remaining : {adata.n_obs}')
    #print(f' {adata.n_obs} observations AFTER Basic Filtering, min_genes={filter_cells_min_genes}')
    sc.pp.filter_cells(adata,min_counts=filter_cells_min_counts)  #  / observations must have min # of coutns
    logger.info(f" Filtering cells sc.pp.filter_cells(adata,min_counts={filter_cells_min_counts})")
    logger.info(f' {adata.n_obs} observations AFTER min_counts Filtering, min_counts={filter_cells_min_counts}')
    #print(f'Filtering cells pp.filter_cells(adata, min_cells=filter_cells_min_counts)  Cells remaining : {adata.n_obs}')
    #print(f'min_counts=filter_cells_min_counts = ',filter_cells_min_counts )
    logger.info(f' {adata.n_vars} features BEFORE Basic Filtering (min_cells and min_counts)')
    #print(f' {adata.n_vars} features BEFORE Basic Filtering (min_cells and min_counts)')
    sc.pp.filter_genes(adata, min_cells=filter_genes_min_cells ) #genes must be present in min # of cells / observations
    logger.info(f" Filtering cells  sc.pp.filter_genes(adata, min_cells={filter_genes_min_cells})")
    logger.info(f' {adata.n_vars} observations AFTER  min_genes Filtering, min_genes={filter_genes_min_cells}')
    #print(f'Filtering genes pp.filter_genes(adata, min_cells=filter_genes_min_cells)  Genes remaining : {adata.n_vars}')
    #print(f' {adata.n_vars} observations AFTER  min_genes Filtering, min_genes={filter_genes_min_cells}')
    sc.pp.filter_genes(adata, min_counts=filter_genes_min_counts ) #genes must have min # of counts for gene to be kept
    logger.info(f" Filtering cells  sc.pp.filter_genes(adata, min_counts={filter_genes_min_counts})")
    logger.info(f' {adata.n_vars} observations AFTER  min_genes Filtering, min_counts={filter_genes_min_counts}')
    #print(f'Filtering genes pp.filter_genes(adata, min_counts=filter_genes_min_counts)  Genes remaining :  {adata.n_vars}')
    #print(f' {adata.n_vars} observations AFTER  min_genes Filtering, min_counts={filter_genes_min_counts}')
    return



def filter_cells_by_anotated_QC_gene(
    adata: anndata.AnnData | None = None,
    filter_ncount: bool | None = True,
    n_genes_bycounts: int = 10000,
    filter_pct_mt: bool | None = True,
    percent_mt: int = 20,
    over_percent_mt: int = 0,
    filter_pct_ribo: bool | None = False,
    percent_ribo: int = 100,
    over_percent_ribo: int = 0,
    filter_pct_hb: bool | None = False,
    percent_hb: int = 100,
    over_percent_hb: int = 0,
    filter_pct_malat1=False,
    percent_malat1=100,
    over_percent_malat1=0,
    **parameters: Any):
    """  Remove cells that have too many mitochondrial genes expressed or too many total counts:
  #### code

                                    
    """
    logger.info(f' {filter_ncount} keep cells with less than {n_genes_bycounts} (n_genes_bycounts) detected genes ')
    logger.info(f' {filter_pct_mt} keep cells with less than {percent_mt} (percent_mt) mitochondiral gene counts ')
    logger.info(f' {filter_pct_mt} keep cells with greater than {over_percent_mt} (percent_mt) mitochondiral gene counts ')
    logger.info(f' {filter_pct_ribo} keep cells with less than {percent_ribo} (percent_ribo) ribosomal protein gene counts ')
    logger.info(f' {filter_pct_ribo} keep cells with greater than {over_percent_ribo} (percent_ribo) ribosomal protein gene counts ')
    logger.info(f' {filter_pct_hb} keep cells with less than {percent_hb} (percent_hb) hemoglobin protein gene counts ')
    logger.info(f' {filter_pct_hb} keep cells with greater than {over_percent_hb} (percent_ribo) ribosomal protein gene counts ')
    logger.info(f' {filter_pct_malat1} keep cells with less than {percent_malat1} (percent_hb) hemoglobin protein gene counts ')
    logger.info(f' {filter_pct_malat1} keep cells with greater than {over_percent_malat1} (percent_ribo) ribosomal protein gene counts ')


    # Actually do the filtering by slicing the `AnnData` object.
    logger.info(f" {adata.n_obs} observations BEFORE QC gene Filtering (n_genes_bycounts,percent_mt, percent_ribo,percent_hb,and percent_malat1)")
    #print(f" {adata.n_obs} observations BEFORE QC gene Filtering (min_genes and min_counts)")
    if filter_ncount ==True:
        adata = adata[adata.obs.n_genes_by_counts <= n_genes_bycounts, :].copy()  # by n_genes_bycounts
        logger.info(f" {adata.n_obs} observations AFTER Filtering min detected genes {n_genes_bycounts} ")
        #print(f" {adata.n_obs} observations AFTER  Filtering min detected genes {n_genes_bycounts} ")
    if filter_pct_mt ==True:
        adata = adata[adata.obs.pct_counts_mt <= percent_mt, :].copy()   # by percent_mt
        logger.info(f" {adata.n_obs} observations AFTER Filtering max percent_mt {percent_mt} ")
        #print(f" {adata.n_obs} observations AFTER  Filtering max percent_mt {percent_mt} ")
        adata = adata[adata.obs.pct_counts_mt >= over_percent_mt, :].copy()    # by percent_mt
        logger.info(f" {adata.n_obs} observations AFTER Filtering min over_percent_mt {over_percent_mt} ")
        #print(f" {adata.n_obs} observations AFTER  Filtering min percent_mt {over_percent_mt} ")
    if filter_pct_ribo ==True:
        adata = adata[adata.obs.pct_counts_ribo <= percent_ribo, :].copy()   # by percent_ribo
        logger.info(f" {adata.n_obs} observations AFTER Filtering min percent_ribo {percent_ribo} ")
        #print(f" {adata.n_obs} observations AFTER Filtering min percent_ribo {percent_ribo} ")
        adata = adata[adata.obs.pct_counts_ribo >= over_percent_ribo, :].copy()    # by percent_ribo
        logger.info(f" {adata.n_obs} observations AFTER Filtering min over_percent_ribo {over_percent_ribo} ")
        #print(f" {adata.n_obs} observations AFTER  Filtering min percent_mt {over_percent_mt} ")
    if filter_pct_hb ==True:
        adata = adata[adata.obs.pct_counts_hb <= percent_hb, :].copy()    # by percent_hb
        logger.info(f" {adata.n_obs} observations AFTER Filtering min percent_hb {percent_hb} ")
        #print(f" {adata.n_obs} observations AFTER Filtering min percent_hb {percent_hb} ")
        adata = adata[adata.obs.pct_counts_hb >= over_percent_hb, :].copy()    # by percent_hb
        logger.info(f" {adata.n_obs} observations AFTER Filtering min over_percent_hb {over_percent_hb} ")
        #print(f" {adata.n_obs} observations AFTER Filtering min over_percent_hb {over_percent_hb} ")
    if filter_pct_malat1 ==True:
        adata = adata[adata.obs.pct_counts_malat1 <= percent_malat1, :].copy()    # by percent_malat1
        logger.info(f" {adata.n_obs} observations AFTER Filtering min percent_malat1 {percent_malat1} ")
        #print(f" {adata.n_obs} observations AFTER Filtering min percent_malat1 {percent_malat1} ")
        adata = adata[adata.obs.pct_counts_malat1 >= over_percent_malat1, :].copy()    # by over_percent_malat1
        logger.info(f" {adata.n_obs} observations AFTER Filtering min over_percent_malat1 {over_percent_malat1} ")
        #print(f" {adata.n_obs} observations AFTER Filtering min over_percent_malat1 {over_percent_malat1} ")
    return adata


def remove_genes(
    adata: anndata.AnnData | None = None,
    remove_MALAT1: bool | None = False,
    remove_MT: bool | None = False,
    remove_HB: bool | None = False,
    remove_RP_SL: bool | None = False,
    remove_MRP_SL: bool | None = False,
    **parameters: Any
):
    """ ################################# Remove Filter out genes with ""techincal bias""
    
     
    """
    ### Remove gene sets  on off switches
    logger.info(f'\n ####################################################  remove_genes')
    logger.info(f'remove_MALAT1 : {remove_MALAT1}')
    logger.info(f'remove_MT : {remove_MT}')
    logger.info(f'remove_HB : {remove_HB}')
    logger.info(f'remove_RP_SL : {remove_RP_SL}')
    logger.info(f'remove_MRP_SL : {remove_MRP_SL} ')
    logger.info(f" {adata.n_vars} features BEFORE Filtering for specific genes and {adata.n_obs} observations") 
    #print(f'####################################################  remove_genes')
    #print(f'remove_MALAT1 : {remove_MALAT1}')
    #print(f'remove_MT : {remove_MT}')
    #print(f'remove_HB : {remove_HB}')
    #print(f'remove_RP_SL : {remove_RP_SL}')
    #print(f'remove_MRP_SL : {remove_MRP_SL} ')
    #print(f" {adata.n_vars} features BEFORE Filtering for specific genes and {adata.n_obs} observations") ")
    nothing = adata.var_names.str.startswith('NO_GENES_HAVE_THIS_NAME')
    remove = np.add(nothing, nothing)
    #print(len((nothing)))
    if remove_MALAT1==True:
        malat1 = adata.var_names.str.startswith('MALAT1')
        remove = np.add(remove, malat1)
    # we need to redefine the mito_genes since they were first 
    # calculated on the full object before removing low expressed genes.
    if remove_MT==True:
        mito_genes = adata.var_names.str.startswith('MT-')
        remove = np.add(remove,mito_genes)
    if remove_HB==True:
        #hb_genes = (adata.var_names.str.startswith('HB')& ~adata.var_names.str.contains(("HBEGF"))  & ~adata.var_names.str.contains(("HBS1L"))  & ~adata.var_names.str.contains(("HBP1"))) # HBEGF,HBS1L, HBP1 not a hemeoglobin genes 
        # "^HB[^(P)" changed to "^HB[^(P)(S)" and  & ~adata_test.var_names.str.contains(("HBEGF")) added to remove HBS1L and HBEGF which are NOT memoglobin genes
        hb_genes = (adata.var_names.str.contains(("^HB[^(P)(S)]")) & ~adata.var_names.str.contains(("HBEGF")) )
        remove = np.add(remove,hb_genes)
    if remove_RP_SL==True:
        RP_SL_genes = adata.var_names.str.startswith(("RPS","RPL"))
        remove = np.add(remove,RP_SL_genes )
    if remove_MRP_SL==True:
        MRP_SL_genes = adata.var_names.str.startswith(("MRPS","MRPL"))
        remove = np.add(remove,MRP_SL_genes )    
    keep = np.invert(remove)
    adata = adata[:,keep].copy()
    logger.info(f" {adata.n_vars} features AFTER Filtering for specific genes and {adata.n_obs} observations") 
    #print(f" {adata.n_vars} features AFTER Filtering for specific genes and {adata.n_obs} observations") 
    return adata    

### filter funcitons END


# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]