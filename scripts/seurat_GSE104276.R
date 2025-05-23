setwd("C:/Abbie/research/seurat/GSE104276/")

library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(Matrix)
library(ggplot2)

## example <- readRDS("GSE104276.CPM_example.rds") # saved seurat object

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
ggsave(filename = "output/scale_and_PCA/GSE104276_VizDimLoadings.jpg", height = 7, width = 12, quality = 50)
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave(filename = "output/scale_and_PCA/GSE104276_DimPlot_PCA.jpg", height = 7, width = 12, quality = 50)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)
ggsave(filename = "output/scale_and_PCA/GSE104276_DimHeatmap.jpg")


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)
ggsave(filename = "output/scale_and_PCA/GSE104276_ElbowPlot.jpg", height = 7, width = 12, quality = 50)

#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = "output/clustering/GSE104276_DimPlot_UMAP.jpg", height = 7, width = 12, quality = 50)
saveRDS(seurat_obj, file = "output/saved_seurat_obj/GSE104276_seurat_obj")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# A tibble: 7,244 × 7
# Groups:   cluster [9]
#         p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
#         <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
#   1 3.04e-104       1.11 0.895 0.592 4.16e-100 0       LDHB     
# 2 1.10e- 81       2.35 0.432 0.111 1.51e- 77 0       CCR7     
# 3 4.20e- 79       1.10 0.848 0.407 5.75e- 75 0       CD3D     
# 4 2.53e- 48       2.12 0.336 0.108 3.46e- 44 0       PRKCQ-AS1
# 5 1.40e- 47       1.21 0.625 0.358 1.91e- 43 0       NOSIP    
# 6 3.99e- 41       1.93 0.316 0.109 5.47e- 37 0       LEF1     
# 7 5.75e- 37       1.37 0.42  0.188 7.89e- 33 0       PIK3IP1  
# 8 1.12e- 32       2.42 0.185 0.045 1.53e- 28 0       FHIT     
# 9 3.67e- 32       1.88 0.259 0.087 5.03e- 28 0       MAL      
# 10 2.99e- 30       2.23 0.155 0.033 4.11e- 26 0       NELL2    
# # ℹ 7,234 more rows
# # ℹ Use `print(n = ...)` to see more rows

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
