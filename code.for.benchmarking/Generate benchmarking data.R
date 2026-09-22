generate_benchamarking_dataset <- function(SampleID, Celltype, PercRare, Ndataset) {
  library(data.table)
  library(Seurat)
  library(dplyr)
  
  load(paste("./data/",SampleID,"/",SampleID,".RData",sep = ""))
  load(paste("./data/",SampleID,"/",SampleID,"_",Celltype,".RData",sep = ""))
  
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
      
      save(P1.seurat.obj, file=paste("data/",SampleID,"/",exprimentID, ".seurat.obj.RData",sep=""))
      
      rm(P1.seurat.obj)
      gc()
    }
  }
}


setwd('/home/jinxiuyuan/Proj_scCellFishing/')
PercRare <- c(0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.0025, 0.00125)
Ndataset <- c(1:10)


for (j in 1:length(Celltype)) {
  generate_benchamarking_dataset(SampleID, Celltype = Celltype[j], PercRare, Ndataset)
}

