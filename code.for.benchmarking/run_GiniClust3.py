#!/usr/bin/env python3

import scanpy as sc
import numpy as np
import anndata
import time
import pandas as pd
import giniclust3 as gc
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
    "GiniClust3"
)

os.makedirs(output_dir, exist_ok=True)


## -----------------------------
## Load expression matrix
## -----------------------------
adataRaw = sc.read_csv(input_file, first_column_names=True)
			
sc.pp.filter_cells(adataRaw,min_genes=3)
sc.pp.filter_genes(adataRaw,min_cells=200)

adataSC=anndata.AnnData(X=adataRaw.X.T,obs=adataRaw.var,var=adataRaw.obs)
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)


## -----------------------------
## Run GiniClust3
## -----------------------------
start_time = time.time()
gc.gini.calGini(adataSC) ###Calculate Gini Index
adataGini=gc.gini.clusterGini(adataSC,neighbors=3) ###Use higher value of neighbor in larger dataset. Recommend (5:15)

gc.fano.calFano(adataSC) ###Calculate Fano factor
adataFano=gc.fano.clusterFano(adataSC) ###Cluster based on Fano factor

consensusCluster={}
consensusCluster['giniCluster']=np.array(adataSC.obs['rare'].values.tolist())
consensusCluster['fanoCluster']=np.array(adataSC.obs['fano'].values.tolist())
gc.consensus.generateMtilde(consensusCluster) ###Generate consensus matrix
gc.consensus.clusterMtilde(consensusCluster) ###Cluster consensus matrix

end_time = time.time()
elapsed_time = end_time - start_time

time_data = {'start_time': [start_time], 'end_time': [end_time], 'elapsed_time': [elapsed_time]}
df = pd.DataFrame(time_data)

time_file = os.path.join(
    output_dir,
    f"time.{experimentID}.csv"
)
df.to_csv(time_file, index=False)

cluster_file = os.path.join(
    output_dir,
    f"{experimentID}.seurat.obj.csv"
)

np.savetxt(
    cluster_file,
    consensusCluster["finalCluster"],
    delimiter="\t",
    fmt="%s"
)



