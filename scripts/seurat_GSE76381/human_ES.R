setwd("C:/Abbie/research/seurat/GSE76381/")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)

#### 1. Setup the Seurat Object ####

raw <- read.table("GSE76381_ESMoleculeCounts.cef.txt", fill = TRUE) 

### 1a. Create metadata ###

metadata <- as.data.frame(t(raw[1:4, ]))
rownames(metadata) <- metadata[, 2]
colnames(metadata) <- metadata[1, ]
metadata <- metadata[-1, ]

### 1b. Create data matrix ###

tmp1 <- raw
colnames(tmp1) <- tmp1[2, ]
rownames(tmp1) <- tmp1[, 1]
tmp1 <- tmp1[-(1:5), -1]
tmp2 <- as.matrix(tmp1) # convert to char matrix
matrix <- matrix(as.numeric(tmp2), 
                 ncol = ncol(tmp2)) # convert to num matrix
colnames(matrix) <- colnames(tmp1)
rownames(matrix) <- rownames(tmp1)

### 1c. Create seurat obj ###

seurat_obj <- CreateSeuratObject(
  counts = as(matrix, "dgCMatrix"),
  meta.data = metadata,
  project = "human_ES"
)
# An object of class Seurat 
# 18539 features across 1715 samples within 1 assay 
# Active assay: RNA (18539 features, 0 variable features)
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
ggsave(filename = "output/scale_and_PCA/human_ES_VizDimLoadings.jpg", height = 7, width = 12, quality = 50)
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave(filename = "output/scale_and_PCA/human_ES_DimPlot_PCA.jpg", height = 7, width = 12, quality = 50)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
ggsave(filename = "output/scale_and_PCA/human_ES_DimHeatmap.jpg")


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)
ggsave(filename = "output/scale_and_PCA/human_ES_ElbowPlot.jpg", height = 7, width = 12, quality = 50)

#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = "output/clustering/human_ES_DimPlot_UMAP.jpg", height = 7, width = 12, quality = 50)
saveRDS(seurat_obj, file = "output/saved_seurat_obj/human_ES_seurat_obj")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# A tibble: 9,652 × 7
# # Groups:   cluster [9]
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene 
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>
#   1 2.18e-145       3.26 0.802 0.137 4.04e-141 0       NFIA 
# 2 5.76e- 88       1.39 0.571 0.098 1.07e- 83 0       TTR  
# 3 2.64e- 87       2.86 0.559 0.094 4.89e- 83 0       HTR2C
# 4 4.49e- 83       3.01 0.58  0.124 8.32e- 79 0       EPHA3
# 5 1.57e- 79       3.17 0.485 0.077 2.91e- 75 0       ZIC4 
# 6 1.75e- 71       1.96 0.846 0.436 3.24e- 67 0       SLIT2
# 7 3.32e- 70       3.14 0.414 0.056 6.15e- 66 0       LRP1B
# 8 4.35e- 70       1.77 0.923 0.8   8.06e- 66 0       SPARC
# 9 5.39e- 70       2.07 0.846 0.482 9.99e- 66 0       WLS  
# 10 1.76e- 68       2.09 0.861 0.517 3.25e- 64 0       FOS  
# # ℹ 9,642 more rows
# # ℹ Use `print(n = ...)` to see more rows

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
