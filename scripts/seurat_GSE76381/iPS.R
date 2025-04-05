setwd("C:/Abbie/research/seurat/GSE76381/")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)

#### 1. Setup the Seurat Object ####

raw <- read.table("GSE76381_iPSMoleculeCounts.cef.txt", fill = TRUE) 

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
  project = "iPS"
)
# An object of class Seurat 
# 14726 features across 337 samples within 1 assay 
# Active assay: RNA (14726 features, 0 variable features)
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
ggsave(filename = "output/scale_and_PCA/iPS_VizDimLoadings.jpg", height = 7, width = 12, quality = 50)
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave(filename = "output/scale_and_PCA/iPS_DimPlot_PCA.jpg", height = 7, width = 12, quality = 50)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)
ggsave(filename = "output/scale_and_PCA/iPS_ElbowPlot.jpg", height = 7, width = 12, quality = 50)

#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = "output/clustering/iPS_DimPlot_UMAP.jpg", height = 7, width = 12, quality = 50)
saveRDS(seurat_obj, file = "output/saved_seurat_obj/iPS_seurat_obj")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# A tibble: 711 × 7
# # Groups:   cluster [3]
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene          
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>         
#   1 6.38e-17       1.05 0.963 0.714  9.40e-13 0       SYT13         
# 2 9.21e-13       4.18 0.451 0.114  1.36e- 8 0       CARTPT        
# 3 1.02e-10       1.31 0.735 0.4    1.50e- 6 0       MIR1244-3-loc4
# 4 1.06e-10       1.92 0.549 0.211  1.57e- 6 0       NR2F2         
# 5 3.12e-10       1.38 0.679 0.354  4.59e- 6 0       MIR1244-3-loc3
# 6 7.75e-10       1.14 0.42  0.109  1.14e- 5 0       ELMO1         
# 7 9.94e-10       1.15 0.747 0.417  1.46e- 5 0       ATP6V0C       
# 8 1.33e- 9       1.42 0.568 0.234  1.97e- 5 0       ONECUT2       
# 9 1.47e- 9       1.08 0.735 0.411  2.16e- 5 0       TSPYL4        
# 10 3.22e- 9       3.09 0.228 0.017  4.73e- 5 0       AMIGO2        
# # ℹ 701 more rows
# # ℹ Use `print(n = ...)` to see more rows

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
