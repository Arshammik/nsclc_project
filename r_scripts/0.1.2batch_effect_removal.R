suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GEOquery))

pbmc <- readRDS("~/Documents/single_cell/nsclc_project/results/pbmc(before_batch_effect_removal).rds")

obj.list <- SplitObject(object = pbmc, split.by = "stage")
obj.list

for (i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <-  FindVariableFeatures(object =  obj.list[[i]])
}

features <- SelectIntegrationFeatures(object.list = obj.list)

# find integrarion anchors (CCA) and (RPCA)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = 'rpca')

anchors <- readRDS("~/Desktop/anchors.rds")
seurat.integrated <- IntegrateData(anchorset = anchors)
saveRDS(seurat.integrated, "~/Documents/single_cell/nsclc_project/results/seurat.integrated.rds")