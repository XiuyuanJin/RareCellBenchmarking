generate_benchamarking_dataset <- function(SampleID, Celltype, PercRare, Ndataset) {
  library(data.table)
  library(Seurat)
  library(dplyr)
  
  load(paste("/data/",SampleID,".RData",sep = ""))
  load(paste("/data/",SampleID,"_",Celltype,".RData",sep = ""))
  
  for (i in PercRare) {
    for (n in Ndataset) {
      
      set.seed(n)
      RareCell  <- colnames(get(paste0(Celltype,".seurat.obj")))
      rmNumRare <- round((dim(seurat.obj)[2]*i-length(RareCell))/(i-1))
      NumRare   <- length(RareCell) - rmNumRare
      rmRare    <- sample(RareCell,rmNumRare)
      
      P1.seurat.obj <- seurat.obj[,setdiff(colnames(seurat.obj),rmRare)]
      exprimentID   <- paste(Celltype,".",i,".",n,sep="")
      print(exprimentID)
      
      save(P1.seurat.obj, file=paste("/results/",exprimentID, ".seurat.obj.RData",sep=""))
      
      rm(P1.seurat.obj)
      gc()
    }
  }
}


PercRare <- 0.05
Ndataset <- 1

SampleID <- "KidneyCellAtlas"
Celltype <- c("ConneTubule")

for (j in 1:length(Celltype)) {
  generate_benchamarking_dataset(SampleID, Celltype = Celltype[j], PercRare, Ndataset)
}

