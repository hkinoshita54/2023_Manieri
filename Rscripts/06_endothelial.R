####################
# description ----
# proceed from seu_03_SCT-harmony_2.1.RDS
# endothelial subset


####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(Seurat)
library(gprofiler2)


####################
# load data, set parameters ----
seu <- readRDS(file = "RDSfiles/seu_03_SCT-harmony_2.1.RDS")

file_no = "06"
description = "ec"
npcs = 15
path = file.path("plots", description, paste0("npc", as.character(npcs)))
dir.create(path = path, recursive = TRUE)
RDSfile = paste0("RDSfiles/", "seu_", file_no, "_", description, "_npc", as.character(npcs), ".RDS")


####################
# subset stromal cells ----
seu <- subset(seu, idents = c("Endothelial"))


####################
# split by orig.ident, lognormalization with harmony integration ----
seu <- DietSeurat(seu)
seu[["SCT"]] <- NULL
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
seu <- SCTransform(seu, variable.features.n = 3000)    # used 3000 var.features
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, 
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, reduction = "harmony", resolution = 0.6, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "tumor") + NoAxes()
ggsave("tumor.png", path = path, width = 5, height = 5, units = "in", dpi = 150)


####################
#### feature plots ----
files <- list.files(path = "gene_set/endothelial_annotation/", full.names = TRUE)
features <- lapply(files, FUN = readLines) %>% unlist()
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

