setwd("/home/jinxiuyuan/Proj_scCellFishing/")

source("software/StemID-master/RaceID2_StemID_class.R")
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
## Run RaceID2
## -----------------------------
start_time  <- Sys.time()
sc          <- SCseq(read.counts)
sc          <- filterdata(sc)
sc          <- clustexp(sc)
sc          <- findoutliers(sc)
final       <- data.frame(CELLID = names(sc@cpart), cluster = sc@cpart)
end_time    <- Sys.time()

compute_time   <- as.numeric(difftime(end_time, start_time, units = "secs")) / 60


## -----------------------------
## Calculate performance metrics
## -----------------------------      
rare.cell   <- colnames(subset(P1.seurat.obj, cell_type %in% clusters[[Celltype]]))
major.cell  <- setdiff(colnames(P1.seurat.obj),rare.cell)

rare.clst    <- names(table(final[final$CELLID %in% rare.cell,]$cluster))[which.max(table(final[final$CELLID %in% rare.cell,]$cluster))]
Raceid.clst  <- subset(final, cluster %in% rare.clst)
Raceid.major <- subset(final, cluster != rare.clst)

TP <- length(intersect(rownames(Raceid.clst),rare.cell))
FP <- length(intersect(rownames(Raceid.clst),major.cell))      
TN <- length(intersect(rownames(Raceid.major),major.cell))
FN <- length(intersect(rownames(Raceid.major),rare.cell))

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
write.csv(result,file = paste0(output_dir, "/RaceID2_", Celltype, ".csv"),row.names = TRUE)

