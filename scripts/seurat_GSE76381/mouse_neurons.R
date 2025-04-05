setwd("C:/Abbie/research/seurat/GSE76381/")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)

#### 1. Setup the Seurat Object ####

raw <- read.table("GSE76381_MouseAdultDAMoleculeCounts.cef.txt", sep = ",", fill = TRUE) 

### 1a. Create metadata ###

metadata <- as.data.frame(t(raw[1:4, ]))
rownames(metadata) <- metadata[, 3]
colnames(metadata) <- c("CEF", "Dataset", "CELL_ID", "Cell_type")
metadata$CEF[2:nrow(metadata)] <- metadata$CEF[1:(nrow(metadata) - 1)] # shift col 'CEF' down by 1
metadata$Dataset[2:nrow(metadata)] <- metadata$Dataset[1:(nrow(metadata) - 1)] # shift col 'Dataset' down by 1
metadata <- metadata[-(1:2), ]

### 1b. Create data matrix ###

tmp1 <- raw
colnames(tmp1) <- tmp1[3, ]
tmp1 <- tmp1[!duplicated(tmp1[, 1]), ] # remove duplicate genes
rownames(tmp1) <- tmp1[, 1]
tmp1 <- tmp1[-(1:4), -(1:2)]
tmp2 <- as.matrix(tmp1) # convert to char matrix
matrix <- matrix(as.numeric(tmp2), 
                 ncol = ncol(tmp2)) # convert to num matrix
colnames(matrix) <- colnames(tmp1)
rownames(matrix) <- rownames(tmp1)

### 1c. Create seurat obj ###

seurat_obj <- CreateSeuratObject(
  counts = as(matrix, "dgCMatrix"),
  meta.data = metadata,
  project = "mouse_neurons"
)
# An object of class Seurat 
# 18220 features across 243 samples within 1 assay 
# Active assay: RNA (18220 features, 0 variable features)
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

# Visualize feature-feature relationships
plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#### 3. Normalizing the data ####

seurat_obj <- NormalizeData(
  seurat_obj, normalization.method = "LogNormalize", 
  scale.factor = 10000
)


#### 4. Identification of highly variable features (feature selection) ####

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#### 5. Scaling the data ####

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)


#### 6. Perform linear dimensional reduction ####

seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

# Examine and visualize PCA results a few different ways
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")
ggsave(filename = "output/scale_and_PCA/mouse_neurons_VizDimLoadings.jpg", height = 7, width = 12, quality = 50)
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave(filename = "output/scale_and_PCA/mouse_neurons_DimPlot_PCA.jpg", height = 7, width = 12, quality = 50)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)
ggsave(filename = "output/scale_and_PCA/mouse_neurons_ElbowPlot.jpg", height = 7, width = 12, quality = 50)

#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = "output/clustering/mouse_neurons_DimPlot_UMAP.jpg", height = 7, width = 12, quality = 50)
saveRDS(seurat_obj, file = "output/saved_seurat_obj/mouse_neurons_seurat_obj")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# # A tibble: 1,840 × 7
# # Groups:   cluster [4]
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
#   1 1.61e-16       1.07 0.988 0.988  2.93e-12 0       Resp18
# 2 4.21e-16       1.06 1     0.969  7.67e-12 0       Gng3  
# 3 1.10e-14       1.07 1     0.975  2.00e-10 0       Tmsb10
# 4 1.67e-12       1.07 0.976 0.963  3.05e- 8 0       Gapdh 
# 5 2.12e-12       1.11 0.951 0.87   3.86e- 8 0       Bex4  
# 6 2.45e-10       1.09 0.951 0.758  4.47e- 6 0       Calb1 
# 7 8.97e-10       1.77 0.622 0.267  1.63e- 5 0       Cthrc1
# 8 1.07e- 9       1.17 0.854 0.727  1.96e- 5 0       Smim18
# 9 1.29e- 9       1.01 0.915 0.807  2.36e- 5 0       Cystm1
# 10 1.81e- 6       2.53 0.341 0.106  3.31e- 2 0       Rpph1 
# # ℹ 1,830 more rows
# # ℹ Use `print(n = ...)` to see more rows

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
