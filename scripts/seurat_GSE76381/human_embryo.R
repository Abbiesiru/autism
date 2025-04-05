setwd("C:/Abbie/research/seurat/GSE76381/")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)

#### 1. Setup the Seurat Object ####

raw <- read.table("GSE76381_EmbryoMoleculeCounts.cef.txt", fill = TRUE) 

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
  project = "human_embryo"
)
# An object of class Seurat 
# 19531 features across 1977 samples within 1 assay 
# Active assay: RNA (19531 features, 0 variable features)
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
ggsave(filename = "output/scale_and_PCA/human_embryo_VizDimLoadings.jpg", height = 7, width = 12, quality = 50)
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
ggsave(filename = "output/scale_and_PCA/human_embryo_DimPlot_PCA.jpg", height = 7, width = 12, quality = 50)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)
ggsave(filename = "output/scale_and_PCA/human_embryo_ElbowPlot.jpg", height = 7, width = 12, quality = 50)

#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
DimPlot(seurat_obj, reduction = "umap")
ggsave(filename = "output/clustering/human_embryo_DimPlot_UMAP.jpg", height = 7, width = 12, quality = 50)
saveRDS(seurat_obj, file = "output/saved_seurat_obj/human_embryo_seurat_obj")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
obj.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
# # A tibble: 13,287 × 7
# # Groups:   cluster [12]
# p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
# <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
#   1 1.19e-149       3.07 0.906 0.346 2.32e-145 0       CELF4    
# 2 1.08e-136       2.57 0.898 0.339 2.10e-132 0       MAPT-loc1
# 3 4.88e-135       2.89 0.822 0.248 9.54e-131 0       HMP19    
# 4 4.44e-122       2.27 0.914 0.421 8.67e-118 0       RTN1     
# 5 4.44e-116       2.12 0.984 0.694 8.67e-112 0       MEG3     
# 6 2.56e-109       2.19 0.836 0.287 5.00e-105 0       INA      
# 7 2.21e-108       2.95 0.674 0.166 4.32e-104 0       MYT1L    
# 8 1.39e- 97       1.55 0.982 0.482 2.72e- 93 0       STMN2    
# 9 1.39e- 95       1.87 0.903 0.533 2.72e- 91 0       DPYSL3   
# 10 4.52e- 94       1.75 0.935 0.661 8.82e- 90 0       NCAM1    
# # ℹ 13,277 more rows
# # ℹ Use `print(n = ...)` to see more rows

obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()
