# Core scverse libraries
from __future__ import annotations
import re

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — renders to file only

import anndata as ad

# Data retrieval
import scanpy as sc
import scanpy.external as sce

# Reproducibility
import scvi
from scvi.external import SysVI
import torch
import scanorama

import celltypist
from celltypist import models


class ScanpyCleaner:
    """
    A pipeline wrapper for standard scanpy preprocessing.
    """
    def __init__(self, adata: sc.AnnData, project_name: str = "Dataset"):
        self.adata = adata
        self.project_name = project_name
        # Keep a copy of raw counts in a layer
        self.adata.layers["counts"] = self.adata.X.copy()

    def run_qc(self, min_genes: int = 500, min_cells: int = 3, min_counts: int = 1000, mito_prefix: str = "MT-", max_mito_pct: float = 20.0):
        """
        Filters cells/genes and calculates mitochondrial metrics.
        """
        print(f"[{self.project_name}] Running QC...")
        
        # Basic filtering
        sc.pp.filter_cells(self.adata, min_genes=min_genes)
        sc.pp.filter_genes(self.adata, min_cells=min_cells)
        sc.pp.filter_cells(self.adata, min_counts=min_counts)

        
        # Calculate QC metrics
        self.adata.var["mt"] = self.adata.var_names.str.startswith(mito_prefix)
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
        
        # Filter by mitochondrial content
        self.adata = self.adata[self.adata.obs["pct_counts_mt"] < max_mito_pct, :].copy()
        print(f"[{self.project_name}] Post-QC: {self.adata.n_obs} cells and {self.adata.n_vars} genes.")

    def normalize(self, target_sum: float = 1e6):
        """
        Normalizes and log-transforms the data.
        """
        print(f"[{self.project_name}] Normalizing and Log-transforming...")
        sc.pp.normalize_total(self.adata, target_sum=target_sum)
        sc.pp.log1p(self.adata)

    def identify_features(self, n_top_genes: int = 2000, flavor: str = "seurat"):
        """
        Identifies Highly Variable Genes (HVGs).
        """
        print(f"[{self.project_name}] Identifying HVGs...")
        sc.pp.highly_variable_genes(self.adata, n_top_genes=n_top_genes, flavor=flavor)
        
    def scale_and_pca(self, max_value: float = 10.0):
        """
        Scales data and runs PCA.
        """
        print(f"[{self.project_name}] Scaling and running PCA...")
        # Store log-normalized data before scaling if needed for visualization
        self.adata.raw = self.adata
        sc.pp.scale(self.adata, max_value=max_value)
        sc.tl.pca(self.adata, svd_solver="arpack")

    def full_pipeline(self):
        """
        Runs the baseline cleaning steps in sequence.
        """
        self.run_qc()
        self.normalize()
        self.identify_features()
        self.scale_and_pca()
        print(f"[{self.project_name}] Pipeline Complete.")
        return self.adata
