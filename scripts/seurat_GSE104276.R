setwd("/Users/abbiesiru/Desktop/research/seurat/GSE104276/")

library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(Matrix)
library(ggplot2)

example <- readRDS("GSE104276.CPM_example.rds") # saved seurat object

#### 1. Setup the Seurat Object ####

tmp1 <- read_excel("GSE104276_all_pfc_2394_UMI_count_NOERCC.xlsx")
tmp1 <- as.data.frame(tmp1)
rownames(tmp1) <- tmp1[, 1]
tmp1 <- tmp1[, -1]
tmp2 <- read.delim("cell_info.GSE104276.txt", header = TRUE, sep = "\t")
matrix <- as.matrix(tmp1[, colnames(tmp1) %in% rownames(tmp2)])
metadata <- tmp2[colnames(matrix), ]

seurat_obj <- CreateSeuratObject(
  counts = as(matrix, "dgCMatrix"),
  meta.data = metadata,
  project = "GSE104276"
)
# An object of class Seurat 
# 24153 features across 2391 samples within 1 assay 
# Active assay: RNA (24153 features, 0 variable features)
# 1 layer present: counts


#### 2. Standard pre-processing workflow ####

### 2a. QC and selecting cells for further analysis ###

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") # none found

# Visualize QC metrics as a violin plot
VlnPlot(
  seurat_obj, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  ncol = 3
)

# Visualize QC metrics as a violin plot
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
