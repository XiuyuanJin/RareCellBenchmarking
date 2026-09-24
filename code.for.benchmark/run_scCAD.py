#!/usr/bin/env python3

from pathlib import Path
import sys
import time
import os
import numpy as np
import pandas as pd


# -----------------------------
# scCAD source
# -----------------------------
software_dir = Path(
    "/code/code.for.benchmarking/Software_source_code"
)

source_file = software_dir / "scCAD.py"

sys.path.insert(0, str(software_dir))
import scCAD


## -----------------------------
## Demo parameters
## -----------------------------
SampleID = "KidneyCellAtlas"
Celltype = "ConneTubule"
PrecRare = 0.05
DatasetID = 1

rare_h = 0.01

## Example experiment ID:
## ConneTubule.0.05.1
experimentID = f"{Celltype}.{PrecRare}.{DatasetID}"
print(experimentID)

input_file = Path("/data") / f"{experimentID}.seurat.obj.csv"

output_dir = "/results"


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
    save_path="/results/"
)

print("scCAD analysis finished.")