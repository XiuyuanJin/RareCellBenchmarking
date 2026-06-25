setwd("/home/jinxiuyuan/Proj_scCellFishing/")

library(CIARA)
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

read.counts <- GetAssayData(P1.seurat.obj,assay = "RNA",slot = "counts") %>% as.matrix()


## -----------------------------
## Run CIARA
## -----------------------------
start_time   <- Sys.time()
data_seurat  <- cluster_analysis_integrate_rare(read.counts, exprimentID, 0.1, 5, 30)
data_cluster <- as.vector(data_seurat$seurat_clusters)

norm_matrix <- GetAssayData(data_seurat,assay = "RNA",slot = "data") %>% as.matrix()
knn_matrix  <- as.matrix(data_seurat@graphs$RNA_nn)

background  <- get_background_full(norm_matrix, threshold = 1, n_cells_low = 3, n_cells_high = 20)

result      <- CIARA(norm_matrix, knn_matrix, background, cores_number = 1, p_value = 0.001, local_region = 1, approximation = FALSE)
ciara_genes <- row.names(result)[result[, 1] < 1]

data_ciara  <- cluster_analysis_integrate_rare(read.counts, exprimentID, 0.01, 5, 30, feature_genes = ciara_genes)
end_time      <- Sys.time()


compute_time   <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60

final_cluster <- merge_cluster(data_cluster, data_ciara$seurat_clusters, max_number = 20)
final_cluster <- as.data.frame(final_cluster)
colnames(final_cluster) <- "cluster"
rownames(final_cluster) <- colnames(P1.seurat.obj)
final_cluster$cluster   <- gsub("-step_2", "", final_cluster$cluster)


## -----------------------------
## Calculate performance metrics
## -----------------------------      
rare.cell  <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell <- setdiff(colnames(P1.seurat.obj),rare.cell)

rare.clst   <- names(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))[which.max(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))]
ciara.clst  <- subset(final_cluster, cluster %in% rare.clst)
ciara.major <- subset(final_cluster, cluster != rare.clst)

TP <- length(intersect(rownames(ciara.clst),rare.cell))
FP <- length(intersect(rownames(ciara.clst),major.cell))      
TN <- length(intersect(rownames(ciara.major),major.cell))
FN <- length(intersect(rownames(ciara.major),rare.cell))

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
write.csv(result,file = paste0(output_dir, "/CIARA_", Celltype, ".csv"),row.names = TRUE)

