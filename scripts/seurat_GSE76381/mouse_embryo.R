setwd("C:/Abbie/research/seurat/GSE76381/")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)

#### 1. Setup the Seurat Object ####

raw <- read.table("GSE76381_MouseEmbryoMoleculeCounts.cef.txt", fill = TRUE) 

### 1a. Create metadata ###

metadata <- as.data.frame(t(raw[1:6, ]))
rownames(metadata) <- metadata[, 3]
colnames(metadata) <- metadata[1, ]
metadata <- metadata[-1, ]

### 1b. Create data matrix ###

tmp1 <- raw
colnames(tmp1) <- tmp1[3, ]
rownames(tmp1) <- tmp1[, 1]
tmp1 <- tmp1[-(1:7), -1]
tmp2 <- as.matrix(tmp1) # convert to char matrix
matrix <- matrix(as.numeric(tmp2), 
                 ncol = ncol(tmp2)) # convert to num matrix
colnames(matrix) <- colnames(tmp1)
rownames(matrix) <- rownames(tmp1)

### 1c. Create seurat obj ###

seurat_obj <- CreateSeuratObject(
  counts = as(matrix, "dgCMatrix"),
  meta.data = metadata,
  project = "mouse_embryo"
)
# An object of class Seurat 
# 24378 features across 1907 samples within 1 assay 
# Active assay: RNA (24378 features, 0 variable features)
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
ggsave(filename = "output/scale_and_PCA/mouse_embryo_VizDimLoadings.jpg", height = 7, width = 12, quality = 50)
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave(filename = "output/scale_and_PCA/mouse_embryo_DimPlot_PCA.jpg", height = 7, width = 12, quality = 50)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)
ggsave(filename = "output/scale_and_PCA/mouse_embryo_ElbowPlot.jpg", height = 7, width = 12, quality = 50)

#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = "output/clustering/mouse_embryo_DimPlot_UMAP.jpg", height = 7, width = 12, quality = 50)
saveRDS(seurat_obj, file = "output/saved_seurat_obj/mouse_embryo_seurat_obj")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# # A tibble: 15,358 × 7
# # Groups:   cluster [11]
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene  
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr> 
#   1 9.61e-117       1.22 0.998 0.896 2.34e-112 0       Tubb2a
# 2 3.70e-107       1.34 0.989 0.624 9.02e-103 0       Snhg11
# 3 3.48e- 96       1.92 0.777 0.252 8.48e- 92 0       Grm5  
# 4 1.80e- 91       1.03 1     0.852 4.39e- 87 0       Stmn2 
# 5 8.66e- 87       1.29 0.954 0.612 2.11e- 82 0       Acot7 
# 6 3.60e- 86       1.54 0.945 0.523 8.78e- 82 0       Snap25
# 7 8.45e- 85       1.74 0.805 0.325 2.06e- 80 0       Fabp3 
# 8 1.28e- 81       1.59 0.867 0.41  3.12e- 77 0       Ndrg4 
# 9 2.34e- 81       1.03 1     0.783 5.71e- 77 0       Gap43 
# 10 5.45e- 81       1.07 0.998 0.744 1.33e- 76 0       Stmn3 
# # ℹ 15,348 more rows
# # ℹ Use `print(n = ...)` to see more rows

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
