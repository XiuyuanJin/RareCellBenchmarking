
source("/code/code.for.benchmarking/Software_source_code/RaceID3_StemID2_class.R")
library(data.table)
library(Seurat)
library(dplyr)


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

read.counts <- as.data.frame(
  as.matrix(
    GetAssayData(
      P1.seurat.obj,
      assay = "RNA",
      slot = "counts"
    )
  )
)


## -----------------------------
## Run RaceID3
## -----------------------------
start_time  <- Sys.time()
sc          <- SCseq(read.counts)
sc          <- filterdata(sc)
sc          <- clustexp(sc)
sc          <- findoutliers(sc)
sc          <- rfcorrect(sc)
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
write.csv(result,file = paste0("/results/RaceID3_", Celltype, ".csv"),row.names = TRUE)

