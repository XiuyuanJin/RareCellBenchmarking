library(Seurat)
library(RSpectra)
library(EDGE)
library(dplyr)
library(FNN)
library(igraph)


## -----------------------------
## Demo parameters
## -----------------------------
SampleID <- "KidneyCellAtlas"
Celltype <- "ConneTubule"
clusters  <- list(ConneTubule = "kidney connecting tubule epithelial cell")

PrecRare <- 0.05
DatasetID <- 1

## Example experiment ID:
## ConneTubule.0.05.1
experimentID <- paste(Celltype, PrecRare, DatasetID, sep = ".")
print(experimentID)


input_file <- paste0("/data/",experimentID, ".seurat.obj.RData")


custom_defs <- EDGE::endr_defs

## -----------------------------
## Load Seurat object
## -----------------------------
load(input_file)
        
dat_log <- t(as.matrix(
  GetAssayData(P1.seurat.obj, assay = "RNA", slot = "data")
))

## -----------------------------
## Run EDGE
## -----------------------------
start_time <- Sys.time()
        
simu.endr <- EDGE::endr(
  as.data.frame(dat_log),
  custom_defs
)
knn_res    <- get.knn(simu.endr, k = 15)
        
edges <- cbind(rep(seq_len(nrow(simu.endr)), each = 15), as.vector(t(knn_res$nn.index)))
g     <- graph_from_edgelist(edges, directed = FALSE)
        
louvain_clusters <- cluster_louvain(g)
cluster_labels   <- membership(louvain_clusters)
        
end_time    <- Sys.time()
compute_time <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60
        
final_cluster    <- data.frame(Cluster = as.numeric(cluster_labels))
colnames(final_cluster) <- "cluster"
rownames(final_cluster) <- colnames(P1.seurat.obj)
        
rare.cell      <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell     <- setdiff(colnames(P1.seurat.obj),rare.cell)
        
rare.clst <- names(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))[which.max(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))]
edge.clst  <- subset(final_cluster, cluster %in% rare.clst)
edge.major <- subset(final_cluster, cluster != rare.clst)
        
TP <- length(intersect(rownames(edge.clst),rare.cell))
FP <- length(intersect(rownames(edge.clst),major.cell))      
TN <- length(intersect(rownames(edge.major),major.cell))
FN <- length(intersect(rownames(edge.major),rare.cell))
        
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
write.csv(result,file = paste0("/results/EDGE_", Celltype, ".csv"),row.names = TRUE)

