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
# load data and check features by Vln ----
seu <- readRDS(file = "RDSfiles/seu_03_SCT-harmony_2.1.RDS")
VlnPlot(seu, 
        features = c("Pdgfra", "Acta2", "Rgs5", "Pecam1"), 
        ncol = 1, pt.size = 0)
VlnPlot(seu, 
        features = c("Ptprc", "S100b", "Epcam", "Mki67"), 
        ncol = 1, pt.size = 0)


####################
# annotation ----
seu <- RenameIdents(
  seu, 
  `0` = "PDGFRAlo", `1` = "Endothelial", `2` = "Endothelial", `3` = "Endothelial", `4` = "PDGFRAlo", 
  `5` = "PDGFRAhi", `6` = "PDGFRAhi", `7` = "Myocyte", `8` = "PDGFRAhi", `9` = "PDGFRAlo", 
  `10` = "Myocyte", `11` = "Immune", `12` ="PDGFRAlo", `13` = "Pericyte", `14` = "Endothelial", 
  `15` = "Immune", `16` = "Immune", `17` = "Endothelial", `18` = "Myocyte", `19` = "Endothelial", `20` = "Immune",
  `21` = "Immune", `22` = "Endothelial", `23` = "Glial", `24` = "Endothelial", `25` = "Pericyte",
  `26` = "Immune", `27` = "Epithelial", `28` = "PDGFRAhi", `29` = "Endothelial"
  )
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, levels = c("PDGFRAhi", "PDGFRAlo", "Myocyte", "Pericyte", "Endothelial", "Immune", "Glial", "Epithelial"))
Idents(seu) <- seu$celltype
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu, group.by = "orig.ident") + NoAxes()
DimPlot(seu, group.by = "position") + NoAxes()

####################
# recheck feature plot ----
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
file.copy(from=plots.png.paths, to="plot/SCT_w_harmony/ver2.1_annotation")

saveRDS(seu, file = "RDSfiles/seu_03_SCT-harmony_2.1.RDS")


####################
# check dot plot ----
features = c("Pdgfra", "Acta2", "Myh11", "Rgs5", "Des", "Pecam1", "Ptprc", "S100b", "Epcam", "Mki67")
DotPlot(seu, features = features) + RotatedAxis()