## module imports
import scanpy as sc



def process2scaledPCA(adata,
                   normalize_total_target_sum=1e4,logarithmize=True,
                   filter_HVG=False,HVG_min_mean=0.0125, HVG_max_mean=3, HVG_min_disp=3,
                   regress_mt=False,regress_ribo=False,regress_malat1=False,regress_hb=False,
                   scale=True,scale_max_std_value=None,PCA=True,
                      cell_cycle_score=True,
                   regress_cell_cycle_score=False,HVG_flavor='seurat',HVG_n_top_genes=1500,organism='human',
                      **parameters):
    '''
    
    #### code
    

                      
    '''
    
    ################################## library-size correct  the data:    
    ################################## Logarithmize    
    adata=norm_log(adata,normalize_total_target_sum,logarithmize, **parameters)
    #################################  HVG selection
    if HVG_flavor=='seurat':
        print(f'############### HVG_flavor=seurat')
        if logarithmize==True:
            adata=HVG_selection_log_norm_seurat(adata,filter_HVG,HVG_min_mean, HVG_max_mean, HVG_min_disp, **parameters)
        else:
            print(f'####################################################  warning data not logerized use alterante HVG selection')
    if HVG_flavor=='seurat_v3':
        print(f'############### HVG_flavor=seurat_v3')
        if logarithmize==True:
            print(f'####################################################  warning data  logerized ....seurat_v3 will use layers["counts"]  ')
        print(f'#################################################### seurat_v3 HVG selection HVG_n_top_genes = {HVG_n_top_genes}   ')
        adata=HVG_selection_log_norm_seurat_v3(adata,filter_HVG,HVG_n_top_genes, **parameters)

    #################################  regress_out_anotated_QC_genes
    regress_out_anotated_QC_genes(adata, regress_mt,regress_ribo, regress_malat1,regress_hb, **parameters)
    if scale==True:
        #print(f'####################################################  Now running sc.pp.scale() Scale the data (each gene to unit variance)')
        #sc.pp.scale(adata, max_value=10)  # Scale and Clip values exceeding standard deviation 10.
        scale_func(adata,scale_max_std_value=scale_max_std_value, **parameters)
    print(f'####################################################  calc_cell_cycle_score')
    print(f'we are calc_cell_cycle_score True/Flase = {cell_cycle_score}')
    if cell_cycle_score ==True:
        print(f'running calc_cell_cycle_score(adata)')
        calc_cell_cycle_score(adata,organism=organism)
    print(f'we are regressing out  cell_cycle_score True/Flase = {regress_cell_cycle_score}')
    if regress_cell_cycle_score ==True:
        regress_cell_cycle_score_func(adata)
    if PCA ==True:
        PCA_func(adata, **parameters)
    else:
        print(f' PCA == False ... skipping PCA ')
    return adata

def scale_func(adata,scale_max_std_value=10, **parameters):
    print(f'####################################################  Now running sc.pp.scale() Scale the data (each gene to unit variance)')
    print(f' scale_max_std_value= {scale_max_std_value}')
    sc.pp.scale(adata, max_value=scale_max_std_value)  # Scale and Clip values exceeding standard deviation 10.
    return 
def PCA_func(adata, **parameters):
    print(f'####################################################  Principal component analysis')
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True)
    return 


def norm_log(adata,normalize_total_target_sum=1e4, logarithmize=True, **parameters):
    '''
    #### code
    def norm_log(adata,normalize_total_target_sum=1e4, logarithmize=True, **parameters):
    
    print(f'####################################################  save raw counts (adata.X) to adata.layers["counts"]')
    adata.layers["counts"] = adata.X.copy()  # preserve counts
    print(f'####################################################  library-size correct  the data  target_sum= {normalize_total_target_sum}')
    ################################## library-size correct  the data:
    sc.pp.normalize_total(adata, target_sum=normalize_total_target_sum)    
    if logarithmize==True:
        print(f'####################################################  Logarithmize  the data')
        print(f'############################# to adata.raw save filtered, normalized and logarithmized gene expression and plot')
        sc.pp.log1p(adata)
    else:
         print(f'############################# to adata.raw save filtered and normalized gene expression and plot')
    adata.raw = adata
    return
    
    '''
    print(f'####################################################  save raw counts (adata.X) to adata.layers["counts"]')
    adata.layers["counts"] = adata.X.copy()  # preserve counts
    print(f'####################################################  library-size correct  the data  target_sum= {normalize_total_target_sum}')
    ################################## library-size correct  the data:
    sc.pp.normalize_total(adata, target_sum=normalize_total_target_sum)
    adata.raw = adata.copy()    #### save normalized counts to adata.raw not logerized
    if logarithmize==True:
        print(f'####################################################  Logarithmize  the data')
        print(f'############################# to adata.raw save filtered, normalized and logarithmized gene expression and plot')
        sc.pp.log1p(adata)
    else:
         print(f'############################# to adata.raw save filtered and normalized gene expression and plot')
    #adata.raw = adata
    return adata


        
def HVG_selection_log_norm_seurat(adata,filter_HVG=False,HVG_min_mean=0.0125, HVG_max_mean=3, HVG_min_disp=3, **parameters):
    '''
    ########################  Identify highly-variable genes from log normed data
    
    #### code
    def HVG_selection_log_norm_seurat(adata,filter_HVG=False,HVG_min_mean=0.0125, HVG_max_mean=3, HVG_min_disp=3, **parameters):
    sc.pp.highly_variable_genes(adata, min_mean=HVG_min_mean, max_mean=HVG_max_mean, min_disp=HVG_min_disp)
    print(f'############################# the number of highly varriable gens are = ',sum(adata.var.highly_variable))
    sc.pl.highly_variable_genes(adata) #### plot HVGs
    if filter_HVG==True:
        adata=HVG_removal(adata)
    else:
        print(f' filter_HVG == False ... all genes will be kept ')
    return adata
    
    '''
    sc.pp.highly_variable_genes(adata, min_mean=HVG_min_mean, max_mean=HVG_max_mean, min_disp=HVG_min_disp)
    print(f'############################# the number of highly varriable gens are = ',sum(adata.var.highly_variable))
    sc.pl.highly_variable_genes(adata) #### plot HVGs
    if filter_HVG==True:
        adata=HVG_removal(adata)
    else:
        print(f' filter_HVG == False ... all genes will be kept ')
    return adata

def HVG_selection_log_norm_seurat_v3(adata,filter_HVG=False,HVG_n_top_genes=1500,**parameters):
    '''
    ########################  Identify highly-variable genes from log normed data
    #### code
    def HVG_selection_log_norm_seurat_v3(adata,filter_HVG=False,HVG_n_top_genes=1500,**parameters):
    sc.pp.highly_variable_genes(adata, n_top_genes=HVG_n_top_genes)
    print(f'############################# the number of highly varriable gens are = ',sum(adata.var.highly_variable))
    sc.pl.highly_variable_genes(adata) #### plot HVGs
    if filter_HVG==True:
        adata=HVG_removal(adata)
    else:
        print(f' filter_HVG == False ... all genes will be kept ')
    return adata
    '''
    sc.pp.highly_variable_genes(adata, n_top_genes=HVG_n_top_genes)
    print(f'############################# the number of highly varriable gens are = ',sum(adata.var.highly_variable))
    sc.pl.highly_variable_genes(adata) #### plot HVGs
    if filter_HVG==True:
        adata=HVG_removal(adata)
    else:
        print(f' filter_HVG == False ... all genes will be kept ')
    return adata

def HVG_removal(adata,filter_HVG=True,**parameters):
    '''
    #### code 
    def HVG_removal(adata):
    print(f' Before  filtering for highly_variable genes : number of Cells {adata.n_obs}, number of genes {adata.n_vars}')
    print(f' filter_HVG = True ... only highly_variable gene will be kept ')
    adata = adata[:, adata.var.highly_variable] # Keep only highly variable genes
    print(f' AFTER  filtering for highly_variable genes: number of Cells {adata.n_obs}, number of genes {adata.n_vars}')
    return adata
    adata = adata[:, adata.var.highly_variable]
    '''
    if filter_HVG == True:
        print(f' Before  filtering for highly_variable genes : number of Cells {adata.n_obs}, number of genes {adata.n_vars}')
        print(f' filter_HVG = True ... only highly_variable gene will be kept ')
        adata = adata[:, adata.var.highly_variable] # Keep only highly variable genes
        print(f' AFTER  filtering for highly_variable genes: number of Cells {adata.n_obs}, number of genes {adata.n_vars}')
    else:
        print(f' filter_HVG = False or not True ... all genes will be kept ')
    return adata

    
def regress_out_anotated_QC_genes(adata,regress_mt=False,regress_ribo=False,regress_malat1=False,regress_hb=False,n_jobs=4, **parameters):
    '''
    #### code
    def regress_out_anotated_QC_genes(adata,regress_mt=False,regress_ribo=False,regress_malat1=False,regress_hb=False,n_jobs=4, **parameters):
    print(f'####################################################  regress_out_anotated_QC_genes ')
    ################################# and Regression 
    print(f'we are regressing out  total_counts {regress_mt}')
    print(f'we are regressing out  pct_counts_mt {regress_mt}')
    print(f'we are regressing out  pct_counts_ribo {regress_ribo}')
    print(f'we are regressing out  pct_counts_malat1 {regress_malat1}')
    print(f'we are regressing out  pct_counts_hb {regress_hb}')
    ################ Do the regression 
    print(f'n_jobs= {n_jobs=}')
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
    '''
    print(f'####################################################  regress_out_anotated_QC_genes ')
    ################################# and Regression 
    print(f'we are regressing out  total_counts {regress_mt}')
    print(f'we are regressing out  pct_counts_mt {regress_mt}')
    print(f'we are regressing out  pct_counts_ribo {regress_ribo}')
    print(f'we are regressing out  pct_counts_malat1 {regress_malat1}')
    print(f'we are regressing out  pct_counts_hb {regress_hb}')
    ################ Do the regression 
    print(f'n_jobs= {n_jobs=}')
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



def calc_cell_cycle_score(adata,organism='human'):
    '''
    ################# calculating cell cycle score 
    data should be scaled first if planing on regressing out cell cycle score 
    
    #### code 
    def calc_cell_cycle_score(adata):
    print(f'############# WARNING data should be scaled first if planing on regressing out cell cycle score')
    # Import cell cycle list and split into s_genes and g2m_genes
    s_genes=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1',         'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    cell_cycle_genes=g2m_genes+s_genes
    print(f' there are {len(s_genes)} s_genes   {len(g2m_genes)} g2m_genes  {len(cell_cycle_genes)} cell_cycle_genes')
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    print(f' there are {len(cell_cycle_genes)} cell_cycle_genes in the dataset')    
     ## do scoring 
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
   # plot the cell cycle scores 
    sc.pl.violin(adata, ['S_score', 'G2M_score'],jitter=0.4,rotation=45)
    return 
    '''

    print(f'############# WARNING data should be scaled first if planing on regressing out cell cycle score')
    if organism=='human': 
        print ('Organism is human, using  human cell cycle genes')
        # Import cell cycle list and split into s_genes and g2m_genes
        #s_genes=['MCM5','PCNA','TYMS','FEN1','MCM2','MCM4','RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7',
        #  'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 
        # 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45',
        #  'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
        #g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2',
        #  'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2',
        #  'AURKB', 'BUB1',         'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1',
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
        print ('Organism is mouse, using  mouse cell cycle genes')
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
    print(f' there are {len(s_genes)} s_genes   {len(g2m_genes)} g2m_genes  {len(cell_cycle_genes)} cell_cycle_genes')
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    print(f' there are {len(cell_cycle_genes)} cell_cycle_genes in the dataset')    
     ## do scoring 
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
   # plot the cell cycle scores 
    sc.pl.violin(adata, ['S_score', 'G2M_score'],jitter=0.4,rotation=45)
    return 

def regress_cell_cycle_score_func(adata):
    '''
    #################  Regress out effects of cell cylce score
    ############# WARNING data should be scaled before scoring if planing on regressing out cell cycle score
    #### CODE
    def regress_cell_cycle_score(adata):
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    #return adata
    return
    '''
    sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
    #return adata
    return



