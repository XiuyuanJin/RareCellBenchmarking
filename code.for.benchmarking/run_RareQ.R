setwd("/home/jinxiuyuan/Proj_scCellFishing/")

library(RareQ)
library(Seurat)
library(dplyr)


## -----------------------------
## Demo parameters
## -----------------------------
SampleID  <- "PBMC3k"
Celltype  <- "Naive"
clusters  <- list(Naive = "Naive CD4 T")

PrecRare <- 0.05
DatasetID <- 1

## Example experiment ID:
## Naive.0.05.1
experimentID <- paste(Celltype, PrecRare, DatasetID, sep = ".")
print(experimentID)


input_file <- paste0("data/", SampleID, "/",experimentID, ".seurat.obj.RData")
output_dir <- paste0("output/", SampleID, "/", Celltype)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


## -----------------------------
## Load Seurat object
## -----------------------------
load(input_file)

P1.seurat.obj <- NormalizeData(P1.seurat.obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA(npcs = 50) %>% 
  RunUMAP(reduction = "pca",dims = 1:50)

P1.seurat.obj <- FindNeighbors(
  object = P1.seurat.obj,
  k.param = 20,
  compute.SNN = FALSE,
  prune.SNN = 0,
  reduction = "pca",
  dims = 1:50,
  force.recalc = FALSE,
  return.neighbor = TRUE
)


## -----------------------------
## Run RareQ
## -----------------------------
start_time    <- Sys.time()
final_cluster <- FindRare(P1.seurat.obj)
end_time      <- Sys.time()

compute_time   <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60

final_cluster  <- final_cluster %>% as.data.frame()
colnames(final_cluster) <- "cluster"
rownames(final_cluster) <- colnames(P1.seurat.obj)


## -----------------------------
## Calculate performance metrics
## -----------------------------      
rare.cell      <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell     <- setdiff(colnames(P1.seurat.obj),rare.cell)

rare.clst <- names(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))[which.max(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))]
RareQ.clst  <- subset(final_cluster, cluster %in% rare.clst)
RareQ.major <- subset(final_cluster, cluster != rare.clst)

TP <- length(intersect(rownames(RareQ.clst),rare.cell))
FP <- length(intersect(rownames(RareQ.clst),major.cell))      
TN <- length(intersect(rownames(RareQ.major),major.cell))
FN <- length(intersect(rownames(RareQ.major),rare.cell))

precision <- TP/(TP+FP)
recall    <- TP/(TP+FN)
F1.score  <- (2*precision*recall)/(precision+recall)

result <- data.frame(
  Experiment = experimentID,
  NumRare = length(rare.cell),
  PrecRare = PrecRare,
  Time = compute_time,
  precision = precision,
  TP = TP,
  FP = FP,
  TN = TN,
  FN = FN,
  recall = recall,
  F1 = F1.score
)

rownames(result) <- experimentID        


## -----------------------------
## Save result
## -----------------------------
write.csv(result,file = paste0(output_dir, "/RareQ_", Celltype, ".csv"),row.names = TRUE)

