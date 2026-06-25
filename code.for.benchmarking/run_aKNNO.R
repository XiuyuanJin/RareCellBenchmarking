
library(aKNNO)
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
print(exprimentID)


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
  RunUMAP(reduction = "pca",dims = 1:30)


## -----------------------------
## Run aKNNO
## -----------------------------
start_time    <- Sys.time()
P1.seurat.obj <- FindNeighbors_aKNNO(P1.seurat.obj,verbose = F)
P1.seurat.obj <- FindClusters(P1.seurat.obj,graph.name="aKNN_O",verbose=F)
end_time      <- Sys.time()
        
comput_time   <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60
        
final_cluster  <- P1.seurat.obj$aKNN_O_res.0.8 %>% as.data.frame()
colnames(final_cluster) <- "cluster"


## -----------------------------
## Calculate performance metrics
## -----------------------------      
rare.cell      <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell     <- setdiff(colnames(P1.seurat.obj),rare.cell)
        
rare.clst      <- names(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))[which.max(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))]
aKNNO.clst     <- subset(final_cluster, cluster %in% rare.clst)
aKNNO.major    <- subset(final_cluster, cluster != rare.clst)
        
TP <- length(intersect(rownames(aKNNO.clst),rare.cell))
FP <- length(intersect(rownames(aKNNO.clst),major.cell))      
TN <- length(intersect(rownames(aKNNO.major),major.cell))
FN <- length(intersect(rownames(aKNNO.major),rare.cell))
        
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
write.csv(result,file = paste0(output_dir, "/aKNNO_", Celltype, "_demo.csv"),row.names = TRUE)

