# DATASET_class.py
"""
DATASET_class ...
"""

from __future__ import annotations
import logging

# -----------------------------------------------------------------------------
# Standard libs
# -----------------------------------------------------------------------------
import logging
from pathlib import Path
from typing import Any, Dict, Optional, List
import os, platform, shlex, subprocess, tarfile,pathlib, re, requests
from urllib.parse import urlparse
from cgi import parse_header   
import gzip, shutil
# -----------------------------------------------------------------------------
# Third‑party libs
# -----------------------------------------------------------------------------
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from sklearn.metrics import silhouette_score, silhouette_samples

# -----------------------------------------------------------------------------
# Local libs
# -----------------------------------------------------------------------------
import single_cell_python_tools as sctl


# -----------------------------------------------------------------------------
# Logging setup
# -----------------------------------------------------------------------------
# set up logging within the module (not the root logger)
import logging
__name__leaf = __name__.split('.')[-1]
logger = logging.getLogger("sctl." + __name__leaf)

#logger = logging.getLogger(__name__)
# log stream handler set up at package import
# logger.setLevel(logging.INFO)
#console_handler = logging.StreamHandler()
#formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#console_handler.setFormatter(formatter)
#if not logger.hasHandlers():
#    logger.addHandler(console_handler)


class DATASET_class:
    """End-to-end single-cell workflow with a central *parameters* dictionary..."""

    # ------------------------------------------------------------------
    # Init & helpers
    # ------------------------------------------------------------------
    def __init__(self, parameters: Dict[str, Any]= {},**kwargs):
        if not isinstance(parameters, dict):
            raise TypeError("parameters must be a dict – got %s" % type(parameters))
        # 1) load defaults ------------------------
        self.parameters: Dict[str, Any] = parameters or {}
        self._set_default_parameters()
        # 2) overwrite with user‑provided dict ----------
        self._apply_parameter_overrides()
        # 3) output directories & Scanpy figure path ------------------------
        self._set_output_directories()
        self._check_if_output_directories_exist(verbose=True)
        if self.parameters.get('make_empty_output_dirs', False):
            logger.warning("'make_empty_output_dirs' found in parameters and set to True \n output directories will be created empty")
            self._make_output_dirs()
        else: 
            logger.warning('(parent) output directory set and prefix set but not made\n add \'make_empty_output_dirs\': True in parameters dictionary to create them at init')
        # 4) set adata if provided
        self.adata: Optional[anndata.AnnData] = None
        # 5) set paths
        self.download_dest_path: str=None
        self.input_file_path: str = parameters.get("input_file_path", "")
        self.path: str = parameters.get("path", "")
        # 6) # Scanpy global settings
        sc.settings.verbosity = 1
        sc.logging.print_header()
        sc.settings.set_figure_params(dpi=80, facecolor="white")
        sc.settings.n_jobs = self.n_jobs
        # 7) set some default class attributes
        self.leiden_clusters_renamed = False
        self.leiden_clustering_done = False
        # 8) make log entries
        logger.info("sctl_DATASET_class initialized with parameters:")
        logger.debug(f" self.parameters  {self.parameters}")
        logger.info(f" self.download_url {self.download_url}")
        logger.info(f" self.download_output_dir {self.download_output_dir}")
        logger.info(f" self.download_output_filename {self.download_output_filename}")
        logger.info(f" self.input_file_path {self.input_file_path}")
        logger.info(f" self.output_prefix {self.output_prefix}")
        logger.info(f" self.output_dir {self.output_dir}")
        
        



    def _set_default_parameters(self) -> None:
        """Populate *all* defaults in a single place."""
        # <<< identical to your old set_default_parameters, but without I/O >>>
        self.defaults: Dict[str, Any] = {
            # meta
            "download_url":'https://www.ncbi.nlm.nih.gov/geo/download/',
            "download_output_dir":None,
            "download_output_filename":None,
            "input_file_path":None,
            'file_format':'h5ad',
            'dataset_prefix_for_10x_triplets':'',
            "output_prefix":'GSE_',
            "output_dir":None, 
            'make_empty_output_dirs': False,  # make empty output dirs if True or if output_dir is not none
            "n_jobs": 4,

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
            "filter_HVG" : False,

            ###less than filter percent 
            #"n_genes_bycounts" : 7000, #less than filter  # now filter_cells_min_genes
            "percent_mt" : 10, #less than filter
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
            "scale_max_std_value":None, # 10 often used
                
            ####regression on off switches
            "regress_mt" : False,
            "regress_ribo" : False,
            "regress_malat1":False,
            "regress_hb" : False,
            "regress_cell_cycle_score" : False,

            ###clustering parameters for clusters
            "number_of_PC" : 30, ### dataset demensionality 
            "number_of_neighbors" : 20,
            "leiden_res" : 1, #leiden clustering resolution


            # UMAP graph parameters
            'umap_marker_gene':True,
            'umap_marker_gene_list': ['IL7R','CD14','LYZ', 'MS4A1','CD8A','GNLY','NKG7','FCGR3A','MS4A7','FCER1A','CST3','PPBP'],

            #cluster naming parameters
            'rename_cluster': True,
                'new_cluster_names' : [  'CD4 T', 'CD14 Monocytes', 'B','CD8 T',   'FCGR3A Monocytes','NK','Dendritic', 'Megakaryocytes'], 
        }
        # install as attributes -------------------------------------------------
        for k, v in self.defaults.items():
            setattr(self, k, v)

    def _apply_parameter_overrides(self):
        """Overlay user‑passed `parameters` dict onto attributes."""
        for k, v in self.parameters.items():
            setattr(self, k, v)

    
    def _merge(self,  overrides: Dict[str, Any]) -> Dict[str, Any]:
        default = self.parameters
        if not isinstance(default, dict):
            raise TypeError(f"parameters must be a dict – got {type(default)}")
        return {**default, **overrides}
    
    def _set_output_directories(self):
        """set path for output directories for tables and figures. and scanpy figure dir"""
        if not self.output_dir or not self.output_prefix:
            #raise ValueError("output_dir and output_prefix must be set")
            # set default output dir to current working directory and prefix to 'dataset_'
            self.output_dir = os.getcwd()
            self.output_prefix = "dataset_"
            logger.warning(f"output dir and and output prefix not set.\n Using defaults dir: {self.output_dir} and dataset specfic file prefix: {self.output_prefix}")
        base = Path(self.output_dir) / self.output_prefix
        self.tables_dir = base / "tables"
        self.figures_dir = base / "figures"
        #self.tables_dir.mkdir(parents=True, exist_ok=True)
        #self.figures_dir.mkdir(parents=True, exist_ok=True)
        # set scanpy settings for figures directory
        sc.settings.figdir = str(self.figures_dir)

    def _check_if_output_directories_exist(self, verbose: bool | None = False) -> None:
        """Check if output directories for tables and figures exist."""
        if self.tables_dir.exists() and self.figures_dir.exists():
            self.output_directories_exist=True
            if verbose:
                logger.warning(f'output directories exist. self.output_directories_exist set to True') 
        else:
            self.output_directories_exist=False
            if verbose:
                logger.warning("One or more output directories do not exist. self.output_directories_exist set to False")

    def _make_output_dirs(self) -> None:
        """Create output directories for tables and figures."""
        if not self.output_dir:
            #raise ValueError("output_dir must be set")
            logger.warning("output_dir not set. setting now")
            self._set_output_directories()
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"made output directories:\n"
                    f"  tables: {self.tables_dir}\n"
                    f"  figures: {self.figures_dir}")

    # ------------------------------------------------------------------
    # I/O – NEW wget/curl implementation
    # ------------------------------------------------------------------

    def download_data(self, **kwargs):
        """
        Download a remote file using *wget* (Linux/macOS) or *curl* (Windows).
        parameters
        ----------
        download_url : str
            URL to download from.
        download_output_dir : str
            Directory to save the downloaded file to.
        """
        # 1) prepare parameters and filenames ------------------------
        kw = self._merge( kwargs)
        url: str = kw["download_url"]  # must be provided via parameters or kwargs
        download_dir = kw.get("download_output_dir",os.getcwd())
        download_filename = kw.get("download_output_filename",None)
        Path(download_dir).mkdir(parents=True, exist_ok=True)
        def _filename_from_head(url: str) -> str | None:
            """Return filename suggested by Content‑Disposition header, or None."""
            try:
                r = requests.head(url, allow_redirects=True, timeout=10)
                cd = r.headers.get("Content-Disposition")
                if cd:
                    _, params = parse_header(cd)
                    return params.get("filename") or params.get("filename*")
            except requests.RequestException:
                logger.warning(f"Failed to parse filename from headers for URL: {url}")
            return None
        # 2) prepare download command ------------------------
        if platform.system().lower().startswith("win"): # Windows ( never tested on windows)
            if download_filename:
                cmd = f"curl -L {shlex.quote(url)} --create-dirs --output {shlex.quote(str(Path(download_dir) / download_filename))}"
            else:
                cmd = f"curl -L {shlex.quote(url)} --create-dirs --output {shlex.quote(str(Path(download_dir) / url.split('/')[-1]))}"
        else:
            if download_filename:
                cmd = f"wget -q --progress=bar:noscroll -O {shlex.quote(os.path.join(download_dir, download_filename))} {shlex.quote(url)}"
            else:
                cmd = f"wget -nv --content-disposition --trust-server-names -P {shlex.quote(download_dir)} {shlex.quote(url)}"
        # 3) Run download command ------------------------
        logger.info(f"Running download command:\n{cmd}")
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        logger.info("Download complete.")
        # 4) resolve the destination file name ------------------------
        if download_filename:
            self.download_dest_path = os.path.join(download_dir, download_filename)
            logger.info(f"Downloaded to {self.download_dest_path}")
        else:
            filename = _filename_from_head(url)
            logger.info(f"filename inferred from headers: {filename}")
            if filename is None:
                # Fallback to the last part of the URL
                filename = urlparse(url).path.split("/")[-1]
                logger.warning(f"Falling back to filename from URL: {filename}")
            self.download_dest_path = os.path.join(download_dir, filename)
            logger.info(f"Downloaded to {self.download_dest_path} (inferred from headers)")
        # set path to the parent directory of the downloaded file
        self.path = str(self.download_dest_path)
        logger.debug(f"Set self.path to {self.path}")
        return self
    
    # ------------------------------------------------------------------
    # unpack_tar
    # ------------------------------------------------------------------
    def unpack_tar(self, **kwargs):
        """Extract a ``.tar``, ``.tar.gz``, or ``.tgz`` archive.

        Parameters (from kwargs or self.parameters)
        -----------------------------------------
        tar_path        : str | None   defaults to ``self.download_dist_path``
        extract_to      : str | None   defaults to os.getcwd()
        set_as_path     : bool         if *True*, `self.path` is set to the first
                                        file extracted (useful for 10x mtx dirs).
        Returns
        -------
        self
        """
        # 1) prepare parameters and filenames
        kw = self._merge(kwargs)
        # Resolve tar path
        tar_path = Path(kw.get("tar_path", self.download_dest_path))
        if not tar_path.exists():
            logger.error(f"File not found: {tar_path}")
            raise FileNotFoundError(tar_path)
        # Prepare extraction root
        extract_root = kw.get("extract_to") or tar_path.parent
        Path(extract_root).mkdir(parents=True, exist_ok=True)
        logger.info(f"Extracting {tar_path} → {extract_root} …")
        # 2) Perform extraction
        mode = "r:*"  # auto‑detect gzip/bzip2/xz/none
        with tarfile.open(tar_path, mode) as tar:
            tar.extractall(path=extract_root)
            members = [Path(extract_root) / m for m in tar.getnames()]
        # 3) Set self.path to the parent dir or first file (if requested)
        set_as_path = kw.get("set_as_path", False)
        if set_as_path and members:
            self.path = str(members[0])
            logger.info(f"Set self.path to first extracted file: {self.path}")
        else:
            # Set self.path to the parent directory that contains a file with 'barcodes' in its name
            barcode_file_path = None
            for file_path in members:
                if "barcodes" in file_path.name.lower():
                    barcode_file_path = file_path
                    break
            if barcode_file_path:
                self.path = str(barcode_file_path.parent)
                logger.info(f"Set self.path to directory containing 'barcodes': {self.path}")
            else:
                self.path = str(tar_path.parent)
                logger.warning("No 'barcodes' file found. Defaulting self.path to extract root.")
                logger.info(f"Set self.path to archive parent: {self.path}")

        logger.info(f"Extracted {len(members)} files to: {extract_root}")
        logger.info(f"Extracted files: {[str(p) for p in members]}")
        #return [str(p) for p in members]
        return self

    def decompress_downloaded_files(self,
        keep_archives: bool = True,
        gunzip_single_files: bool = True,
        unpack_archive: bool = False,
         **kwargs ):
        """
        Decompress any archives sitting in ``self.parameters["download_output_dir"]``.
        Parameters
        ----------
        keep_archives : bool, default True
            If ``False`` the original archive (``*.tar.gz``, ``*.zip``…) is deleted
            after successful extraction.
        gunzip_single_files : bool, default False
            If ``True`` also expand lone ``*.gz`` files (e.g. ``matrix.mtx.gz``) in
            place; otherwise they are skipped.
        Returns
        -------
        self
            Enables call-chaining.
        """
        out_dir = Path(self.parameters["download_output_dir"]).expanduser().resolve()
        if not out_dir.exists():
            raise FileNotFoundError(f"Download directory {out_dir} does not exist")
        for path in out_dir.iterdir():
            if path.suffixes[-2:] == [".tar", ".gz"] or path.suffix == ".tgz":
                # ---------- tar.gz / tgz ----------
                logger.info(f"Extracting {path.name} …" )
                with tarfile.open(path, "r:gz") as tar:
                    tar.extractall(path=out_dir)
                if not keep_archives:
                    path.unlink()
            elif path.suffix == ".tar" and not unpack_archive:
                # ---------- plain tar ----------
                logger.info(f"Extracting {path.name} …")
                with tarfile.open(path, "r") as tar:
                    tar.extractall(path=out_dir)
                if not keep_archives:
                    path.unlink()
            elif path.suffix == ".zip":
                # ---------- zip ----------
                logger.info(f"Extracting {path.name} …")
                shutil.unpack_archive(path, out_dir, "zip")
                if not keep_archives:
                    path.unlink()
            elif path.suffix == ".gz" and gunzip_single_files:
                # ---------- lone .gz (e.g. matrix.mtx.gz) ----------
                target = path.with_suffix("")  # strip .gz
                logger.info(f"Decompressing {path.name} > {target.name}")
                with gzip.open(path, "rb") as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                if not keep_archives:
                    path.unlink()
            else:
                logger.debug(f"No decompression needed for {path.name}")
        return self
    
    def fixtypo_in_downloaded_file_name(self,
                                                  downloaded_file='',
                                                renamed_downloaded_file_name=''):
        '''Fix the typo in downloaded file name '''
        # import os and shutil to rename the file
        import os
        import shutil
        downloaded_file_dir=self.parameters['download_output_dir']
        # check if the file exists
        if not os.path.exists(os.path.join(downloaded_file_dir, downloaded_file)):
            #raise FileNotFoundError(f"File {downloaded_file} does not exist in {downloaded_file_dir}")
            logger.warning(f"File {downloaded_file} does not exist in {downloaded_file_dir}")
        # check if the new file name exists
        if os.path.exists(os.path.join(downloaded_file_dir, renamed_downloaded_file_name)):
            #raise FileExistsError(f"File {renamed_downloaded_file_name} already exists in {downloaded_file_dir}")
            logger.warning(f"File {renamed_downloaded_file_name} already exists in {downloaded_file_dir}")
        # rename 'downloaded_file' to 'renamed_downloaded_file_name'
        shutil.move(os.path.join(downloaded_file_dir, downloaded_file),
                    os.path.join(downloaded_file_dir,renamed_downloaded_file_name))
        logger.info(f"Renamed {downloaded_file} to {renamed_downloaded_file_name} in {downloaded_file_dir}")
        return self

    # ------------------------------------------------------------------
    # Core pipeline steps 
    # QC
    # Transform
    # Clustering
    # Ploting processed data
    # ------------------------------------------------------------------

    
    def load_data(self, **kwargs):
        """
        Load data from a file.
        Parameters
        ----------
        path : str
            Path to the file to load.
        file_format : str
            Format of the file to load. Can be one of 'h5ad' or "prefix_10x", '10x'.
        dataset_prefix_for_10x_triplets : str
            Prefix for 10X triplets dataset. Used only if file_format is 'prefix_10x'. do not include the underscore at the end.
            if None or '', the prefix will be inferred from the barcodes file in the directory.
        kwargs : dict
            Additional keyword arguments to pass to the loading function.
        Returns
        -------
        """
        kw = self._merge( kwargs)
        file_format = kw.get("file_format", "h5ad")
        dataset_prefix_for_10x_triplets = kw.get("dataset_prefix_for_10x_triplets", None)

        #.) check if outputdirectories exist if not make them 
        self._check_if_output_directories_exist()
        if not self.output_directories_exist:
            self._make_output_dirs()

        if self.path== "" and self.input_file_path == "":
            raise ValueError("path and input_file_path must be set")
        if self.path== "":
            # set self.path to self.input_file_path
            self.path = self.input_file_path            
        if file_format == "h5ad":
            self.adata = sc.read_h5ad(self.path,)
            logger.info(f"Loaded h5ad file from {self.path}")
            ########## set self.loaded_file_format to h5ad
            self.loaded_file_format = "h5ad"
        elif file_format == "prefix_10x":
            if not dataset_prefix_for_10x_triplets:
                # parse the 10X file prefix from the path directory
                files_in_self_path = os.listdir(self.path)
                barcodes_files = [f for f in files_in_self_path if 'barcodes' in f]
                if len(barcodes_files) == 0:
                    raise ValueError(f"No barcodes files found in {self.path}")
                barcode_file = barcodes_files[0]
                self.file_prefix_10x = barcode_file.split('_')[0]+'_'
                logger.debug(f"10X barcode file: {barcode_file} with prefix {self.file_prefix_10x}")
            else:
                self.file_prefix_10x = dataset_prefix_for_10x_triplets+'_'
                logger.debug(f"Using provided 10X prefix: {self.file_prefix_10x}")
            ########## load the 10X data with the prefix
            self.adata = sc.read_10x_mtx(self.path,prefix=self.file_prefix_10x)
            ########## set self.loaded_file_format to prefix_10x
            self.loaded_file_format = "prefix_10x"
            logger.info(f"Loaded 10X data from {self.path} with prefix {self.file_prefix_10x}")
        elif file_format == "10x":
            self.adata = sc.read_10x_mtx(self.path,prefix=None)
            logger.info(f"Loaded 10X data from {self.path} with no prefix")
            ########## set self.loaded_file_format to 10x
            self.loaded_file_format = "10x"
        else:
            raise ValueError(f"Unsupported format: {file_format}")
        self.adata.uns["parameters"] = self.parameters
        # #) make log entries
        logger.info(f"Loaded data with {self.adata.n_obs} observations and {self.adata.n_vars} variables")
        logger.info(f"self.adata \n{self.adata}")
        return self
    
    # ------------------------------------------------------------------
    # QC : Core pipeline step
    # ------------------------------------------------------------------

    def basic_filitering(self, **kwargs):
        """ Basic Filtering """
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.basic_filitering(self.adata,**kw)
        return self


    
    ### annotate_QC_genes and calculate_qc_metrics methods
    def annotate_QC_genes(self , **kwargs):
        """ 
        sctl.pp.annotate_QC_genes(self.adata,organism = organism,)
        """
        # set up the parameters for the function
        kw = self._merge( kwargs)
        organism = kw.get("organism", "human")
        logger.info(f"Annotating QC genes for organism: {organism}")
        # cal the annotate_QC_genes function
        sctl.pp.annotate_QC_genes(self.adata,**kw)
        return self

    
    def calculate_qc_metrics(self, **kwargs):
        """ calculate_qc_metrics"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.calculate_qc_metrics(self.adata)
        return self

    
    def annotate_n_view_adata_raw_counts(self,):
        """  Annotate technical gene groups  and calculate qc metrics"""
        self.annotate_QC_genes()
        self.calculate_qc_metrics()
        self.plot_qc_metrics()
        return self
    
    def plot_qc_metrics(self):
        """ plot_qc_metrics of Annotated technical gene groups  and top 20 highly expressed"""
        sctl.pp.plot_QC_metrics_violin(self.adata)  
        sctl.pp.plot_QC_metrics_scatter(self.adata) 
        #sc.pl.highest_expr_genes(self.adata, n_top=20, )
    
    ### filter cells and remove genes methods
    def filter_cells_by_anotated_QC_gene(self, **kwargs):
        """  Remove cells that have too many mitochondrial genes expressed or too many total counts:""" 
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.filter_cells_by_anotated_QC_gene(self.adata,**kw)
        return self
    def remove_genes(self, **kwargs):
        """ ################################# Remove Filter out genes with ""techincal bias""
        ### Remove gene sets  on off switches
        """ 
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.remove_genes(self.adata,**kw)
        return self
    
    # ------------------------------------------------------------------
    # transform and normalize  : Core pipeline step
    # ------------------------------------------------------------------
    ### transform and normalize methods
    def norm_log(self, **kwargs):
        """  norm_log"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.norm_log(self.adata,**kw)
        return self
    def HVG_selection_log_norm_seurat(self, **kwargs):
        """  HVG_selection_log_norm_seurat"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.HVG_selection_log_norm_seurat(self.adata,**kw)
        return self
    def HVG_selection_log_norm_seurat_v3(self, **kwargs):
        """  HVG_selection_log_norm_seurat_v3"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.HVG_selection_log_norm_seurat_v3(self.adata,**kw)
        return self
    def HVG_removal(self, **kwargs):
        """  HVG_removal"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.HVG_removal(self.adata,**kw)
        return self
    def regress_out_anotated_QC_genes(self, **kwargs):
        """  regress_out_anotated_QC_genes"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.regress_out_anotated_QC_genes(self.adata,**kw)
        return self
    def scale_func(self, **kwargs):
        """  scale_function"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.scale_func(self.adata,**kw)
        return self
    def PCA_func(self, **kwargs):
        """  PCA_func"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.PCA_func(self.adata,**kw)
        return self
    def calc_cell_cycle_score(self, **kwargs):
        """  calc_cell_cycle_score_func"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.calc_cell_cycle_score(self.adata,**kw)
        return self
    def regress_cell_cycle_score_func(self, **kwargs):
        """  regress_out_cell_cycle_score"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.regress_cell_cycle_score_func(self.adata,**kw)
        return self
    # all transform and normalize methods ran together
    def process2scaledPCA(self, **kwargs):
        """  process2scaledPCA"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pp.process2scaledPCA(self.adata,**kw)
        return self
    # ------------------------------------------------------------------
    # clustering   : Core pipeline step
    # ------------------------------------------------------------------
    ### _clustering methods
    def leiden_clustering(self, **kwargs):
        """  leiden_clustering"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.leiden_clustering(self.adata,**kw)
        # set self.leiden_clustering_done to True
        self.leiden_clustering_done = True
        return self
    def rename_leiden_clusters(self, **kwargs):
        """  rename_leiden_clusters"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.rename_leiden_clusters(self.adata,**kw)
        # set self.leiden_clusters_renamed to True
        self.leiden_clusters_renamed = True
        return self
    def leiden_cluster_sil_score(self, **kwargs):
        """  leiden_cluster_sil_score"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.leiden_cluster_sil_score(self.adata,**kw)
        return self
    def silhouette_walk_Largest_drop(self, **kwargs):
        """  silhouette_walk_Largest_drop"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.silhouette_walk_Largest_drop(self.adata,**kw)
        return self
    def silhouette_walk_4_Largest_drops(self, **kwargs):
        """  silhouette_walk_4_Largest_drops"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pp.silhouette_walk_4_Largest_drops(self.adata,**kw)
        return self
    # ------------------------------------------------------------------
    # Ploting processed data   : Core pipeline step
    # ------------------------------------------------------------------
    def marker_gene_umaps(self, **kwargs):
        """  marker_gene_umaps"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        # Custom colormap where zero values are represented by grey
        import matplotlib.colors as mcolors
        cmap = mcolors.ListedColormap(['gray'] + list(plt.cm.viridis(np.linspace(0, 1, 256))))
        # extract parameters from kw
        vmin = kw.get("vmin", 0)
        vmax = kw.get("vmax", 'p98')
        ncols = kw.get("ncols", 6)
        palette= kw.get("palette", sc.pl.palettes.godsnot_102[1:])
        wspace= kw.get("wspace", None)
        title= kw.get("title", None)
        rename_clusters = kw.get("rename_cluster",False)
        # plot the marker genes
        marker_genes = kw.get("marker_genes", self.adata.uns["parameters"]['umap_marker_gene_list'])
        # check if the marker_genes are in the adata.var_names
        if not all(gene in self.adata.var_names for gene in marker_genes):
            missing_genes = [gene for gene in marker_genes if gene not in self.adata.var_names]
            logger.warning(f"Marker genes {missing_genes} are not in the adata.var_names")
            # remove the missing genes from the marker_genes list
            marker_genes = [gene for gene in marker_genes if gene in self.adata.var_names]
            logger.info(f"Using marker genes {marker_genes} for plotting")
        if self.leiden_clusters_renamed:
            additonal_plots=kw.get("additonal_plots", ['leiden', 'Cell_Clusters_Named'])  
        elif self.leiden_clustering_done:
            additonal_plots=kw.get("additonal_plots", ['leiden'])
        else:
            additonal_plots=kw.get("additonal_plots", [])

        sc.pl.umap(self.adata, color=marker_genes + additonal_plots,ncols=ncols,
                   wspace=wspace,title=title,vmax=vmax,vmin=vmin, palette=palette,cmap=cmap)
    
    ### Ploting processed data mehtods
    def silhouette_score_n_plot(self, **kwargs):
        """  silhouette_score_n_plot"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        sctl.pl.silhouette_score_n_plot(self.adata,**kw)
        return self
    def silhouette_score_of_obs_key_n_plot(self, **kwargs):
        """  silhouette_score_of_obs_key_n_plot"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pl.silhouette_score_of_obs_key_n_plot(self.adata,**kw)
        return self
    def plot_batch_obs_key_of_obs_key2(self, **kwargs):
        """  plot_batch_obs_key_of_obs_key2"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pl.plot_batch_obs_key_of_obs_key2(self.adata,**kw)
        return self
    def plot_percent_obs_key2_per_batch_obs_key( **kwargs):
        """  plot_percent_obs_key2_per_batch_obs_key"""
        # set up the parameters for the function
        kw = self._merge( kwargs)
        self.adata=sctl.pl.plot_percent_obs_key2_per_batch_obs_key(self.adata,**kw)
        return self
    


 
 # ------------------------------------------------------------------
# Auto-export: collect every function or class defined *in this file*
# whose name does NOT start with an underscore
# ------------------------------------------------------------------
__all__ = [name for name in dir() if not name.startswith("_")]