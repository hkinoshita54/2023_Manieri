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
# check genes of interest from Dr. Hayakawa ----
features = readLines("gene_set/ENS_related.csv")
for (i in 1:length(features)){
  p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
  print(p)
}
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="plot/SCT_w_harmony/ver2.1_ens_related")
rm(plots.dir.path, plots.png.paths)

DotPlot(seu, features = features) + RotatedAxis()


####################
# check other genes ----
features = readLines("gene_set/BMP_related.txt")
DotPlot(seu, features = features) + RotatedAxis()

features = readLines("gene_set/WNT_related.txt")
DotPlot(seu, features = features) + RotatedAxis()

features = c(readLines("gene_set/FGF_related.txt"),readLines("gene_set/EGF_related.txt"))
DotPlot(seu, features = features) + RotatedAxis()
