####################
# load packages ----
library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(glmGamPoi)
library(SeuratWrappers)
library(SeuratDisk)
library(reticulate)
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1e10)


####################
# load data and convert to Assay 5 ----
seu <- readRDS(file = "RDSfiles/seu_01_DF_0.1.RDS")
seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay5")
Idents(seu) <- "orig.ident"
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
seu <- subset(seu, subset = nCount_RNA > 1000 & nFeature_RNA < 3000 & percent.mt < 12)     # nCount>1000, mt<12 in the paper (the upper limit of nFeature was not clear)

####################
# cluster with integration ----
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

seu <- SCTransform(seu, vars.to.regress = c("percent.mt"), 
                   variable.features.n = 3000    # 2000 in the paper
)

seu <- RunPCA(seu, npcs = 30, verbose = FALSE)    # npcs = 15 in the paper

seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony",
  verbose = FALSE
)

seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, reduction = "harmony", resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, group.by = "orig.ident") + NoAxes()
DimPlot(seu, group.by = "position") + NoAxes()

FeaturePlot(seu,features = "Pdgfra", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Acta2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Myh11", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Rgs5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Des", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Pecam1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Ptprc", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "S100b", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Epcam", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "Mki67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="plot/SCT_w_harmony/ver2.1")

saveRDS(seu, file = "RDSfiles/seu_03_SCT-harmony_2.1.RDS")