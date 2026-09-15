#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import os
sys.path.append('/home/jinxiuyuan/Proj_scCellFishing/software/scCAD')
import scCAD

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
    "scCAD"
)

os.makedirs(output_dir, exist_ok=True)


## -----------------------------
## Load expression matrix
## -----------------------------
data_df = pd.read_csv(input_file, index_col=0)

## The CSV file is expected to be genes x cells.
geneNames = np.array(data_df.index)
cellNames = np.array(data_df.columns)

## scCAD expects cells x genes.
data = data_df.T
data = np.array(data, dtype=float)


## -----------------------------
## Run scCAD
## -----------------------------
result, score, sub_clusters, degs_list = scCAD.scCAD(
    data=data,
    dataName=f"{Celltype}_{PrecRare}_{DatasetID}",
    cellNames=cellNames,
    geneNames=geneNames,
    rare_h=0.01,
    save_path=str(output_dir) + "/"
)

print("scCAD analysis finished.")