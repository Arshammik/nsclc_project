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

mtx2 <- readRDS("~/Documents/single_cell/nsclc_project/data/GSE131907_sparse_matrix.rds")
dim(mtx2)

colsums_mtx1 <- Matrix::colSums(mtx1)
colsums_mtx2 <- Matrix::colSums((mtx2))

mtx1_log2 <- Matrix::Matrix(log2(mtx1 + 1), sparse = TRUE)
colsums_mtx1_log2 <- Matrix::colSums(mtx1_log2)

hist(colsums_mtx1, main = "Original mtx1")
hist(log2(colSums(mtx1) + 1), main = "Log2 Scaled colsums(mtx1)")
hist(colsums_mtx2, main = "Original mtx2")
hist(colsums_mtx1_log2, main = "Log2 Scaled mtx1")

features_mtx1 <- rownames(mtx1)
features_mtx2 <- rownames(mtx2)
features <- intersect(features_mtx1, features_mtx2)

mtx1_log2 <- mtx1_log2[features,]
mtx2 <- mtx2[features,]
dim(mtx1) ; dim(mtx2)

mtx <- cbind(mtx1_log2, mtx2)
saveRDS(mtx, '~/Documents/single_cell/nsclc_project/data/mtx.rds')