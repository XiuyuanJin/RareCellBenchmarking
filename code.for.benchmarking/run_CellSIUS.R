setwd("/home/jinxiuyuan/Proj_scCellFishing/")

library(CellSIUS)
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
  RunPCA() %>% 
  FindNeighbors() %>%
  FindClusters()

norm.counts    <- P1.seurat.obj@assays$RNA@data %>% as.matrix()
cluster        <- P1.seurat.obj$seurat_clusters %>% as.character()
names(cluster) <- colnames(P1.seurat.obj)


## -----------------------------
## Run CellSIUS
## -----------------------------
start_time     <- Sys.time()
CellSIUS.out   <- CellSIUS(mat.norm = norm.counts, 
                           group_id = cluster, 
                           min_n_cells = 10, 
                           min_fc = 2,
                           corr_cutoff = NULL, 
                           iter = 0, 
                           max_perc_cells = 50,
                           fc_between_cutoff = 1, 
                           mcl_path = "/home/jinxiuyuan/mcl/local/bin/mcl")
end_time       <- Sys.time()


comput_time   <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60

Result_List    <- CellSIUS_GetResults(CellSIUS.out = CellSIUS.out)
final_cluster  <- CellSIUS_final_cluster_assignment(CellSIUS.out, group_id = cluster, min_n_genes = 3) %>% as.data.frame()
colnames(final_cluster) <- "cluster"
rownames(final_cluster) <- colnames(P1.seurat.obj)


## -----------------------------
## Calculate performance metrics
## -----------------------------      
rare.cell      <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell     <- setdiff(colnames(P1.seurat.obj),rare.cell)

rare.clst      <- names(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))[which.max(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))]
CellSius.clst  <- subset(final_cluster, cluster %in% rare.clst)
CellSius.major <- subset(final_cluster, cluster != rare.clst)

TP <- length(intersect(rownames(CellSius.clst),rare.cell))
FP <- length(intersect(rownames(CellSius.clst),major.cell))      
TN <- length(intersect(rownames(CellSius.major),major.cell))
FN <- length(intersect(rownames(CellSius.major),rare.cell))

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
write.csv(result,file = paste0(output_dir, "/CellSIUS_", Celltype, ".csv"),row.names = TRUE)

