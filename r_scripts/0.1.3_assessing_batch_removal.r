suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(BatchQC))

merged_Seurat_object <- "~/Documents/single_cell/nsclc_project/results/seurat.integrated.rds"
pbmc_before_batch_effect_removal <- "~/Documents/single_cell/nsclc_project/results/pbmc(before_batch_effect_removal).rds"
merged_pbmc <- readRDS(merged_Seurat_object)
pbmc <- readRDS(pbmc_before_batch_effect_removal)


#scaling and dimensional reductions
merged_pbmc <- ScaleData(object = merged_pbmc)
merged_pbmc <- RunPCA(object = merged_pbmc)
merged_pbmc <- RunUMAP(object = merged_pbmc, dims = 1:50)

png("../results/plots/1.1_DimPlot_pbmc.png", width = 1000, height = 1000)
DimPlot(object = pbmc, reduction = "umap",group.by = "stage", raster = FALSE)
dev.off()

png("../results/plots/1.2_DimPlot_mergedpbmc.png", width = 1000, height = 1000)
DimPlot(object = merged_pbmc, reduction = "umap",group.by = "stage", raster = FALSE)
dev.off()

