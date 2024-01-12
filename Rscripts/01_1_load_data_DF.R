####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
library(DoubletFinder)
options(Seurat.object.assay.version = "v3")

####################
# load data ----
sra_table <- read.delim(file = "data/SraRunTable.txt", sep = ",") # read SraRunTable from GSE224737
ids <- sra_table$Sample.Name[3:6] %>% sort() %>% as.character() # get character vector of sample ids, only for the scRNA-seq of stroma
pos <- factor(c("Antrum", "Antrum", "Corpus", "Corpus")) # to add position to the meta data

seu_list <- list()
for (i in 1:length(ids)){
  mtx <- dir(path = paste0("./data/GSE224737_RAW/", ids[i]), pattern = "matrix", full.names = TRUE)
  features <- dir(path = paste0("./data/GSE224737_RAW/", ids[i]), pattern = "features", full.names = TRUE)
  barcodes <- dir(path = paste0("./data/GSE224737_RAW/", ids[i]), pattern = "barcodes", full.names = TRUE)
  cts <- ReadMtx(mtx = mtx, features = features, cells = barcodes)
  seu <- CreateSeuratObject(counts = cts, project = ids[i], min.cells = 3, min.features = 200)
  seu$position <- pos[i]
  seu_list <- append(seu_list, seu)
}

####################
# filter by percent.mt ----
for (i in seq_along(seu_list)) {
  seu_list[[i]]$percent.mt <- PercentageFeatureSet(seu_list[[i]], pattern = "^mt-")
  seu_list[[i]] <- subset(seu_list[[i]], subset = percent.mt < 12)
}

####################
# DoubletFinder ----
for (i in 1:length(seu_list)) {
  # Pre-process seurat object with standard seurat workflow
  seu <- NormalizeData(seu_list[[i]])
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:10)
  seu <- FindNeighbors(object = seu, dims = 1:10)              
  seu <- FindClusters(object = seu, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(seu, PCs = 1:10)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(0.08 * nrow(seu@meta.data)) ## Assuming 8% doublet formation rate - tailor for your dataset, see below
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  seu <- doubletFinder_v3(seu = seu, 
                          PCs = 1:10, 
                          pK = optimal.pk,
                          nExp = nExp.poi.adj)
  metadata <- seu@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder"
  seu@meta.data <- metadata 
  
  # subset and save
  singlets <- subset(seu, doublet_finder == "Singlet")
  seu_list[[i]] <- singlets
  remove(singlets)
}

####################
# merge the list and save ----
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(ids)], add.cell.ids = ids)
saveRDS(seu, file = "RDSfiles/seu_01_DF.RDS")



