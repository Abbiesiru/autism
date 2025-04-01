setwd("C:/Abbie/research/seurat/GSE104276/")

library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(Matrix)
library(ggplot2)


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
DimPlot(seurat_obj, reduction = "pca") + NoLegend()
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(seurat_obj)


#### 8. Cluster the cells ####

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_obj), 5)


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat_obj, reduction = "umap")
saveRDS(seurat_obj, file = "output/GSE104276_preclustering")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find all markers of cluster 2
cluster2.markers <- FindMarkers(seurat_obj, ident.1 = 2)
head(cluster2.markers, n = 5)
#             p_val avg_log2FC pct.1 pct.2    p_val_adj
# LTB  1.612441e-82  1.2420031 0.983 0.642 2.211301e-78
# IL32 1.159398e-77  1.1329390 0.944 0.472 1.589998e-73
# LDHB 3.451452e-74  1.1199331 0.972 0.607 4.733321e-70
# CD3D 1.368338e-65  0.9953969 0.914 0.438 1.876539e-61
# IL7R 5.971448e-58  1.3300602 0.732 0.333 8.189244e-54

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(seurat_obj, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
#                       p_val avg_log2FC pct.1 pct.2     p_val_adj
# FCGR3A        1.448731e-208   6.893999 0.975 0.037 1.986790e-204
# CFD           2.362385e-199   6.301821 0.937 0.033 3.239775e-195
# IFITM3        4.090907e-199   6.242196 0.975 0.045 5.610270e-195
# RP11-290F20.3 1.025050e-189   6.265284 0.849 0.017 1.405753e-185
# CD68          4.386481e-188   5.549969 0.906 0.033 6.015620e-184

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

cluster0.markers <- FindMarkers(seurat_obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(seurat_obj, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(seurat_obj, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(seurat_obj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
obj.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(seurat_obj, features = top10$gene) + NoLegend()


#### 11. Assigning cell type identity to clusters ####

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "output/images/GSE104276_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(seurat_obj, file = "output/GSE104276_final.rds")
