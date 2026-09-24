library(dplyr)
library(Seurat)
library(SCISSORS)
library(future)


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


## -----------------------------
## Load Seurat object
## -----------------------------
load(input_file)


## -----------------------------
## Run SCISSORS
## -----------------------------
start_time    <- Sys.time()

future::plan(future::sequential)

P1.seurat.obj <- SCISSORS::PrepareData(
  seurat.object = P1.seurat.obj,
  n.HVG = 4000,
  n.PC = 15,
  which.dim.reduc = "umap",
  use.parallel = FALSE,
  random.seed = 629
)
sil_score_df <- ComputeSilhouetteScores(P1.seurat.obj, avg = FALSE)
        
cluster_avg <- sil_score_df %>% group_by(Cluster) %>% summarise(AvgScore = mean(Score))
overall_avg <- mean(cluster_avg$AvgScore)
which.clust <- cluster_avg %>% filter(AvgScore < overall_avg) %>% pull(Cluster) %>% as.character() %>% as.numeric()         
        
plan(sequential)
t_reclust <- ReclusterCells(P1.seurat.obj, 
                            which.clust = which.clust, 
                            merge.clusters = TRUE, 
                            k.vals = c(30, 40, 50), 
                            resolution.vals = c(.2, .3, .4), 
                            n.HVG = 4000, 
                            n.PC = 15, 
                            redo.embedding = TRUE, 
                            use.parallel = FALSE,
                            random.seed = 312)
        
P1.seurat.obj <- IntegrateSubclusters(P1.seurat.obj, reclust.results = t_reclust)

end_time       <- Sys.time()
compute_time    <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60
        
final_cluster  <- FetchData(P1.seurat.obj,  vars = "seurat_clusters") %>% as.data.frame()
colnames(final_cluster) <- "cluster"
rownames(final_cluster) <- colnames(P1.seurat.obj)

## -----------------------------
## Calculate performance metrics
## -----------------------------      
rare.cell      <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell     <- setdiff(colnames(P1.seurat.obj),rare.cell)

rare.clst <- names(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))[which.max(table(final_cluster[rownames(final_cluster) %in% rare.cell, ]))]
SCISSORS.clst  <- subset(final_cluster, cluster %in% rare.clst)
SCISSORS.major <- subset(final_cluster, cluster != rare.clst)

TP <- length(intersect(rownames(SCISSORS.clst),rare.cell))
FP <- length(intersect(rownames(SCISSORS.clst),major.cell))      
TN <- length(intersect(rownames(SCISSORS.major),major.cell))
FN <- length(intersect(rownames(SCISSORS.major),rare.cell))

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
write.csv(result,file = paste0("/results/SCISSORS_", Celltype, ".csv"),row.names = TRUE)

