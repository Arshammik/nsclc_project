library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(ggplot2)
library(Matrix)
library(bigmemory)
library(tictoc)
library(GEOquery)
library(data.table)

#gathering the file paths
file_paths <- paste0("~/Documents/single_cell/nsclc_project/data/GSE189357_RAW/", 
                    list.files(path = "~/Documents/single_cell/nsclc_project/data/GSE189357_RAW/"))

for(i in 1:length(file_paths)){
  name <- sub("(^.*)(...$)", "\\2", file_paths[i])
  data <- Read10X(file_paths[i])
  colnames(data) <- paste0(colnames(data), "_", name)
  assign(name, data)
}

mtx1 <- cbind(TD1, TD2, TD3, TD4, TD5, TD6, TD7, TD8, TD9)
dim(mtx1)

mtx2 <- readRDS("~/Documents/Single Cell/Own Project/Data/GSE131907_Lung_Cancer_normalized_log2TPM_matrix_filtered.rds")
mtx2 <- bigmemory::as.matrix(mtx2)
mtx2 <- Matrix(mtx2, sparse = TRUE)
dim(mtx2)


features_mtx1 <- rownames(mtx1)
features_mtx2 <- rownames(mtx2)
features <- intersect(features_mtx1, features_mtx2)

mtx1 <- mtx1[features,]
mtx2 <- mtx2[features,]
dim(mtx1) ; dim(mtx2)

mtx <- cbind(mtx1, mtx2)
saveRDS(mtx, '~/Documents/Single Cell/Own Project/Data/mtx(GSE189357 & GSE131907).rds')

