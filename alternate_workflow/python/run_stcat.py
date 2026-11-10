import sys
import multiprocessing

class _SyncProcess:
    """Drop-in replacement for multiprocessing.Process that runs synchronously."""
    def __init__(self, target=None, args=()):
        self._target = target
        self._args = args
    def start(self):
        if self._target:
            self._target(*self._args)
    def join(self):
        pass

multiprocessing.Process = _SyncProcess

import os
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ann
import STCAT
import warnings
warnings.filterwarnings('ignore')

input_file = sys.argv[1]
output_file = sys.argv[2]
input_dir = os.path.dirname(input_file)

print(f"Reading matrix from: {input_file}")
adata = sc.read_h5ad(input_file)
print(f"Matrix loaded: {adata.n_obs} cells x {adata.n_vars} genes")

obs_csv = os.path.join(input_dir, "obs_metadata.csv")
if os.path.exists(obs_csv):
    print("Attaching obs metadata...")
    obs = pd.read_csv(obs_csv, index_col=0)
    obs = obs.reindex(adata.obs_names)
    adata.obs = obs

emb_csv = os.path.join(input_dir, "embedding.csv")
if os.path.exists(emb_csv):
    print("Attaching Harmony/PCA embedding for neighbour graph...")
    emb_df = pd.read_csv(emb_csv, index_col=0)
    emb_df = emb_df.reindex(adata.obs_names)
    adata.obsm['X_emb'] = emb_df.values.astype(np.float32)
    print(f"Embedding attached: {adata.obsm['X_emb'].shape}")

print(f"AnnData ready: {adata.n_obs} cells x {adata.n_vars} genes")

import scipy.sparse as sp
mat = sp.csr_matrix(adata.X, dtype=np.float32)
print(f" data dtype: {mat.data.dtype}")
print(f" indices dtype: {mat.indices.dtype}")
print(f" indptr dtype: {mat.indptr.dtype}")

mat.indices = mat.indices.astype(np.int32, copy=False)
mat.indptr = mat.indptr.astype(np.int32, copy=False)
adata.X = mat
del mat
print(f" Matrix normalised to CSR float32 / int32 indices")

if 'neighbors' not in adata.uns:
    if 'X_emb' in adata.obsm:

        print("Pre-computing neighbourhood graph from Seurat embedding (Harmony/PCA)...")
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_emb')
        print(f"Neighbourhood graph ready: {adata.n_obs} cells")
    else:

        print("No embedding found - pre-computing neighbours via sparse PCA (arpack, n_comps=30)...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=30, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
        print(f"Neighbourhood graph ready: {adata.n_obs} cells")
else:
    print("Existing neighbourhood graph detected - skipping pre-computation")

results = STCAT.STCAT(adata)

print("Extracting predictions and freeing memory...")
obs_data = results.obs.copy()
del results, adata
import gc
gc.collect()

for _col in obs_data.columns:
    if not pd.api.types.is_numeric_dtype(obs_data[_col]):
        obs_data[_col] = obs_data[_col].astype(str)

print(f"Writing results (metadata only) to: {output_file}")
minimal = ann.AnnData(
    X = np.zeros((len(obs_data), 1), dtype=np.float32),
    obs = obs_data,
)
minimal.write_h5ad(output_file)
del minimal, obs_data
print("STCAT completed successfully!")
