library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GEOquery)

mtx <- readRDS("~/Documents/single_cell/nsclc_project/data/mtx.rds")
pbmc <- CreateSeuratObject(counts = mtx, project = 'nsclc', min.cells = 2, min.features = 200)
pbmc[["mt.percent"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")

VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "mt.percent"), ncol = 3, raster = FALSE)
pbmc <- subset(pbmc, subset = nFeature_RNA < 5000 & nCount_RNA < 10000 & mt.percent < 5)
rownames_meat <- rownames(pbmc@meta.data)
pbmc@meta.data$barcode <- rownames_meat
pbmc@meta.data$title <- sub("(^.{17})(.*$)","\\2",rownames_meat)


# we have two GEO data sets 
uniqie_rownames <- unique(metadata_rownames)
#GSE131907
GSE131907 <- getGEO(GEO = "GSE131907", GSEMatrix = TRUE)
GSE131907_pheno <- pData(object = GSE131907[[1]])
GSE131907_pheno <- GSE131907_pheno [,c(1,2,8,11,47,48)]
colnames(GSE131907_pheno) <- c("title", "geo_accession", "tissue_source", "tumor_stage", "tissue_origin", "stage")
GSE131907_pheno <- GSE131907_pheno[GSE131907_pheno$title %in% uniqie_rownames,]
GSE131907_pheno[1:11,6] <- "N"

#GSE189357
GSE189357 <- getGEO(GEO = "GSE189357", GSEMatrix = TRUE)
GSE189357_pheno <- pData(GSE189357[[1]])
GSE189357_pheno <- GSE189357_pheno[,c(1,2,10,11,12)]
GSE189357_pheno$title <- sub(" scRNA-seq","",GSE189357_pheno$title)
GSE189357_pheno$stage <- c("IV", "IV", "I", "I", "0", "I", "0", "0", "IV")

sample_stages <- rbind(GSE131907_pheno[,c(1,6)], GSE189357_pheno[,c(1,6)])
sample_stages[37:45,1] <- paste0("1_", sample_stages[37:45,1])

pbmc@meta.data <- merge(pbmc@meta.data, sample_stages, by = "title")
rownames(pbmc@meta.data) <- pbmc@meta.data$barcode

pbmc <- NormalizeData(object = pbmc, normalization.method = "RC", scale.factor = 1e6)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", loess.span = 0.3, mean.function = 2000)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
ElbowPlot(object = pbmc)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- RunUMAP(object = pbmc, dims = 1:20)

DimPlot(object = pbmc, reduction = "umap", group.by = "stage", raster = FALSE)

