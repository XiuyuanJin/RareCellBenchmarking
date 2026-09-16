#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy import sparse
import scanpy as sc
import anndata as ad
from shannonca.dimred import reduce
import time
import os

PROJECT_DIR = "/home/jinxiuyuan/Proj_scCellFishing"

## -----------------------------
## Demo parameters
## -----------------------------
SampleID = "PBMC3k"

Celltype = "Naive"
PrecRare = 0.05
DatasetID = 1

## Example experiment ID:
## Naive.0.05.1
experimentID = f"{Celltype}.{PrecRare}.{DatasetID}"
print(experimentID)

input_file = os.path.join(
    PROJECT_DIR,
    "data",
    SampleID,
    f"{experimentID}.seurat.obj.csv"
)

output_dir = os.path.join(
    PROJECT_DIR,
    "output",
    SampleID,
    "SCA"
)

os.makedirs(output_dir, exist_ok=True)


## -----------------------------
## Load expression matrix
## -----------------------------
df = pd.read_csv(input_file, index_col=0)
df = df.transpose()
df = df.apply(pd.to_numeric, errors="raise")

X = sparse.csr_matrix(df.values.astype(float))
cell_names = df.index.astype(str)
gene_names = df.columns.astype(str)


## -----------------------------
## Run scCAD
## -----------------------------
start_time = time.time()
res = reduce(
    X,
    n_comps=50,       
    iters=1,
    nbhd_size=15,
    metric='euclidean',
    model='wilcoxon',
    chunk_size=1000,
    n_tests='auto',
    keep_scores=True,
    keep_loadings=True,
    keep_all_iters=False
)
			
X_sca = res['reduction'] 

adata = ad.AnnData(X_sca)
adata.obs_names = cell_names

sc.pp.neighbors(
    adata,
    n_neighbors=15,
    use_rep='X', 
    metric='euclidean'
)

sc.tl.leiden(adata, resolution=1.0, key_added='SCA_clusters')

end_time = time.time()
elapsed_time = end_time - start_time

time_file = os.path.join(
    output_dir,
    f"time.{experimentID}.csv"
)

pd.DataFrame({
    "start_time": [start_time],
    "end_time": [end_time],
    "elapsed_time": [elapsed_time]
}).to_csv(time_file, index=False)

cluster_file = os.path.join(
    output_dir,
    f"{experimentID}.SCA.cluster.csv"
)

final_clusters = adata.obs['SCA_clusters']
final_clusters.to_csv(cluster_file)
