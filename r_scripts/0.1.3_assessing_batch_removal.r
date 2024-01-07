suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(BatchQC))

merged_Seurat_object <- "~/Documents/single_cell/nsclc_project/results/seurat.integrated.rds"
merged_pbmc <- readRDS(merged_Seurat_object)
