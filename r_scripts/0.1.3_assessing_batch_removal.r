suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GEOquery))
set.seed(42)

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

pbmc.markers <- FindAllMarkers(merged_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

refrence_20 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

refrence_100 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)

write_csv(x = refrence_20, file = "~/Documents/single_cell/nsclc_project/results/refrence_20.csv")
write_csv(x = refrence_100, file = "~/Documents/single_cell/nsclc_project/results/refrence_100.csv")
saveRDS(merged_pbmc, "~/Documents/single_cell/nsclc_project/results/seurat.integrated(find_markers).rds")