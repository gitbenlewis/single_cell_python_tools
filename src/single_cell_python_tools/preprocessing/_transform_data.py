## module imports
import scanpy as sc
import anndata
from typing import Any, Dict, Optional, List, Union

# set up logging within the module (not the root logger)
import logging
__name__leaf = __name__.split('.')[-1]
logger = logging.getLogger("sctl.pp." + __name__leaf)

def process2scaledPCA(
    adata: anndata.AnnData | None = None,
    normalize_total_target_sum: int = 1e4,
    logarithmize: bool = True,
    filter_HVG: bool = False,
    HVG_min_mean: float = 0.0125,
    HVG_max_mean: float = 3,
    HVG_min_disp: float = 3,
    regress_mt: bool = False,
    regress_ribo: bool = False,
    regress_malat1: bool = False,
    regress_hb: bool = False,
    scale: bool = True,
    scale_max_std_value=None,
    PCA: bool = True,
    cell_cycle_score: bool = True,
    regress_cell_cycle_score: bool = False,
    HVG_flavor: str = 'seurat',
    HVG_n_top_genes: int = 1500,
    organism: str = 'human',
    **parameters: Any):
    '''
    This function performs a series of preprocessing steps on the AnnData object, including normalization, log transformation, highly variable gene selection, regression of unwanted sources of variation, scaling, cell cycle scoring, and PCA.
    '''
    logger.info(f" running process2scaledPCA adata at start \n {adata}")
    ################################## library-size correct the data:
    logger.info(f"process2scaledPCA()-Step 1) library-size correct and Logarithmize (optional) the data ")
    logger.info(f"{normalize_total_target_sum} : normalize_total_target_sum ")
    logger.info(f"{logarithmize} : Logarithmize (optional) the data ")
    #adata=norm_log(adata,normalize_total_target_sum,logarithmize, **parameters)
    adata = norm_log(
    adata,
    normalize_total_target_sum=normalize_total_target_sum,
    logarithmize=logarithmize,
    **parameters)
    #################################  HVG selection
    logger.info(f"Select and annotate highly-variable genes (HVGs) ")
    logger.info(f"HVG selection flavor : HVG_flavor= {HVG_flavor} ")
    if HVG_flavor=='seurat':
        if logarithmize==True:
            #logger.info(f' logarithmize==True')
            adata=HVG_selection_log_norm_seurat(adata,filter_HVG,HVG_min_mean, HVG_max_mean, HVG_min_disp, **parameters)
        else:
            logger.warning(f' warning data not logerized use alterante HVG selection')
    if HVG_flavor=='seurat_v3':
        if logarithmize==True:
            logger.warning(f' warning data is logerized ....seurat_v3 will use layers["counts"]  ')
        logger.info(f'{HVG_n_top_genes} : HVG_n_top_genes, HVG_flavor==seurat_v3 HVG selection ')
        adata=HVG_selection_log_norm_seurat_v3(adata,filter_HVG,HVG_n_top_genes, **parameters)

    #################################  regress_out_anotated_QC_genes
    logger.info(f"process2scaledPCA()-Step 2) regress_out_anotated_QC_genes (optional)")
    regress_out_anotated_QC_genes(adata, regress_mt,regress_ribo, regress_malat1,regress_hb, **parameters)
    logger.info(f"process2scaledPCA()-Step 3) scale the data (optional) : {scale} ")
    if scale==True:
        scale_func(adata,scale_max_std_value=scale_max_std_value, **parameters)
    logger.info(f"process2scaledPCA()-Step 4) cell cycle score (optional)  ")
    logger.info(f" {cell_cycle_score} : calculate cell cycle score (optional) ")
    logger.info(f" {regress_cell_cycle_score} : regressing out cell cycle score (optional) ")
    if cell_cycle_score ==True:
        logger.info(f"running sctl.pp.calc_cell_cycle_score(adata,organism={organism}) ")
        calc_cell_cycle_score(adata,organism=organism)
    if regress_cell_cycle_score ==True:
        regress_cell_cycle_score_func(adata)
    logger.info(f"process2scaledPCA()-Step 5) PCA analysis   ")
    if PCA ==True:
        PCA_func(adata, **parameters)
    else:
        logger.warning(f"PCA == False ... skipping PCA, PCA required for UMAP and tSNE")
    logger.info(f"Done : process2scaledPCA adata AFTER \n {adata}")
    return adata

def scale_func(adata,scale_max_std_value=None, **parameters):
    logger.info(f' Scale the data (each gene to unit variance)')
    logger.info(f"running sc.pp.scale(adata,scale_max_std_value={scale_max_std_value}) ")
    sc.pp.scale(adata, max_value=scale_max_std_value)  # Scale and Clip values exceeding standard deviation 10.
    return 
def PCA_func(adata, **parameters):
    logger.info(f' Principal component analysis, results stored in adata.obsm["X_pca"]')
    logger.info(f"running sc.tl.pca(adata, svd_solver='arpack')")
    logger.info(f"plotting variance ratio sc.pl.pca_variance_ratio(adata, log=True)")
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    return 


def norm_log(
    adata: anndata.AnnData | None = None,
    normalize_total_target_sum: float = 1e4,
    save_counts_layer: bool = True,
    logarithmize: bool = True,
    use_lognorm_for_raw: bool = False,
    **parameters: Any):
    '''
    '''
    logger.info(f'sctl.pp.norm_log() ')
    logger.info(f'norm_log()-Step 1) preserve counts (optional) : {save_counts_layer} ')
    if save_counts_layer==True:
        logger.info(f'save raw counts (adata.X) to adata.layers["counts"]')
        adata.layers["counts"] = adata.X.copy()  # preserve counts
    ################################## library-size correct  the data:
    logger.info(f'norm_log()-Step 2) library-size correct each observation target_sum= {normalize_total_target_sum} ')
    sc.pp.normalize_total(adata, target_sum=normalize_total_target_sum)
    ################################## logarithmize the data (optional) and save to adata.raw:
    if use_lognorm_for_raw==False:
        logger.info(f' using normalized data for adata.raw set becasue use_lognorm_for_raw = {use_lognorm_for_raw}')
        logger.info(f'norm_log()-Step 3) save normalized counts to adata.raw not logerized and not scaled \n adata.raw used for plotting')
        adata.raw = adata.copy()    #### save normalized counts to adata.raw not logerized
        if logarithmize==True:
            logger.info(f'norm_log()-Step 4) logarithmize the data (optional): {logarithmize}')
            sc.pp.log1p(adata)
        else:
            logger.info(f'data not logerized because logarithmize (switch)= {logarithmize}')
    elif use_lognorm_for_raw==True:
        if logarithmize==True:
            logger.info(f' using log(normalized) data for adata.raw set becasue use_lognorm_for_raw = {use_lognorm_for_raw}')
            logger.info(f'norm_log()-Step 3) logarithmize the data (optional): {logarithmize}')
            sc.pp.log1p(adata)
            logger.info(f'norm_log()-Step 4) save lognorm for adata.raw , saving lognormed data (logCP10k) to adata.raw ')
            adata.raw = adata.copy()
        else:
            logger.warning(f'  use_lognorm_for_raw = {use_lognorm_for_raw} and logarithmize = {logarithmize} ... this is not recommended, adata.raw will be saved as logCP10k but not logerized')
            logger.info(f'norm_log()-Step 3) save normalized counts to adata.raw not logerized and not scaled \n adata.raw used for plotting')
            adata.raw = adata.copy()
    return adata



        
def HVG_selection_log_norm_seurat(
    adata: anndata.AnnData | None = None,
    filter_HVG: bool = False,
    HVG_min_mean: float = 0.0125,
    HVG_max_mean: float = 3,
    HVG_min_disp: float = 3,
    **parameters: Any):
    '''
    Select highly variable genes
    '''
    sc.pp.highly_variable_genes(
    adata,
    min_mean=HVG_min_mean,
    max_mean=HVG_max_mean,
    min_disp=HVG_min_disp)
    n_HVGs = sum(adata.var.highly_variable)
    logger.info(f'{n_HVGs} : number of highly varriable genes ')
    sc.pl.highly_variable_genes(adata) #### plot HVGs
    if filter_HVG==True:
        adata=HVG_removal(adata)
    else:
        logger.info(f' filter_HVG == {filter_HVG} ... all genes will be kept ')
    return adata

def HVG_selection_log_norm_seurat_v3(
    adata: anndata.AnnData | None = None,
    filter_HVG: bool = False,
    HVG_n_top_genes: int = 1500,
    **parameters: Any
):
    '''
    Select highly variable genes
    '''
    sc.pp.highly_variable_genes(adata, n_top_genes=HVG_n_top_genes)
    n_HVGs = sum(adata.var.highly_variable)
    logger.info(f'{n_HVGs} : number of highly varriable genes ')
    sc.pl.highly_variable_genes(adata) #### plot HVGs
    if filter_HVG==True:
        adata=HVG_removal(adata)
    else:
        logger.info(f' filter_HVG == {filter_HVG}  ... all genes will be kept ')
    return adata

def HVG_removal(
    adata: anndata.AnnData | None = None,
    filter_HVG: bool = True,
    **parameters: Any
):
    '''
    Remove lowly variable genes from the dataset.
    '''
    n_HVGs = sum(adata.var.highly_variable)
    logger.info(f'{n_HVGs} : number of highly varriable genes ')
    if filter_HVG == True:
        logger.info(f' {adata.n_vars} features BEFORE filtering for highly_variable genes')
        logger.info(f' filter_HVG = {filter_HVG} ... only highly_variable gene will be kept ')
        adata = adata[:, adata.var.highly_variable] # Keep only highly variable genes
        logger.info(f' {adata.n_vars} features AFTER filtering for highly_variable genes')
    else:
        logger.info(f' filter_HVG = {filter_HVG} or not True ... all genes will be kept ')
    return adata

    
def regress_out_anotated_QC_genes(
    adata: anndata.AnnData | None = None,
    regress_mt: bool = False,
    regress_ribo: bool = False,
    regress_malat1: bool = False,
    regress_hb: bool = False,
    n_jobs: int = 1,
    **parameters: Any
):
    '''
    Regress out annotated QC genes
    '''
    logger.info(f'####################################################  regress_out_anotated_QC_genes ')
    ################################# and Regression 
    logger.info(f'{regress_mt} (optional) regressing out  total_counts ')
    logger.info(f'{regress_mt} (optional) regressing out  pct_counts_mt')
    logger.info(f'{regress_ribo} (optional) regressing out  pct_counts_ribo')
    logger.info(f'{regress_malat1} (optional) regressing out  pct_counts_malat1')
    logger.info(f'{regress_hb} (optional) regressing out  pct_counts_hb ')
    ################ Do the regression 
    logger.info(f'n_jobs= {n_jobs=}')
    if regress_mt ==True:
        # by total_counts
        sc.pp.regress_out(adata, ['total_counts' ],n_jobs=n_jobs)
    if regress_mt ==True:
        # by percent_mt
        sc.pp.regress_out(adata, ['pct_counts_mt' ],n_jobs=n_jobs)
    if regress_ribo ==True:
        # by percent_ribo
        sc.pp.regress_out(adata, ['pct_counts_ribo' ],n_jobs=n_jobs)
    if regress_malat1 ==True:
        # by percent_hb
        sc.pp.regress_out(adata, ['pct_counts_malat1'],n_jobs=n_jobs)
    if regress_hb ==True:
        # by percent_hb
        sc.pp.regress_out(adata, ['pct_counts_hb' ],n_jobs=n_jobs)
    #return adata
    return



def calc_cell_cycle_score(
    adata: anndata.AnnData | None = None,
    organism: str = 'human',
    **parameters: Any
):
    '''
    calculate cell cycle score based on a s_genes and g2m_genes list
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    '''
    logger.info(f'############# WARNING data should be scaled first if planing on regressing out cell cycle score')
    if organism=='human': 
        logger.info('Organism is human, using  human cell cycle genes')
        # Import cell cycle list and split into s_genes and g2m_genes
        #s_genes=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7',
        #  'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 
        # 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45',
        #  'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
        #g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2',
        #  'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2',
        #  'AURKB', 'BUB1',    'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1',
        #  'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8',
        #  'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
        s_genes  = [
        "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6","CDCA7",
        "DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN",
        "WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45",
        "CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B",
        "BRIP1","E2F8"
        ]
        g2m_genes = [
            "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2","NUF2",
            "CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2","CKAP2L","CKAP2",
            "AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3","HN1",
            "CDC20","TTK","CDC25C","KIF2C","RANGAP1","NCAPD2","DLGAP5","CDCA2","CDCA8",
            "ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP5","CENPE","CTCF",
            "NEK2","G2E3","GAS2L3","CBX5","CENPA"
        ]
    elif organism == "mouse":
        # Same orthologs, mouse symbols
        logger.info('Organism is mouse, using  mouse cell cycle genes')
        s_genes  = [
            "Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7",
            "Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp","Rad51ap1","Gmnn",
            "Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2","Rad51","Rrm2","Cdc45",
            "Cdc6","Exo1","Tipin","Dscc1","Blm","Casp8ap2","Usp1","Clspn","Pola1","Chaf1b",
            "Brip1","E2f8"
        ]
        g2m_genes = [
            "Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2","Nuf2",
            "Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2","Ckap2l","Ckap2",
            "Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp","Cdca3","Hn1",
            "Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5","Cdca2","Cdca8",
            "Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5","Cenpe","Ctcf",
            "Nek2","G2e3","Gas2l3","Cbx5","Cenpa"
        ]
    cell_cycle_genes=g2m_genes+s_genes
    logger.info(f' there are {len(s_genes)} s_genes   {len(g2m_genes)} g2m_genes  {len(cell_cycle_genes)} cell_cycle_genes')
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    logger.info(f' there are {len(cell_cycle_genes)} cell_cycle_genes in the dataset')    
    ## do scoring 
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    # plot the cell cycle scores 
    sc.pl.violin(adata, ['S_score', 'G2M_score'],jitter=0.4,rotation=45)
    return 

def regress_cell_cycle_score_func(
    adata: anndata.AnnData | None = None,
    **parameters: Any
):
    '''
    Regress out cell cycle scores from the data
    '''
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    return





# ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]