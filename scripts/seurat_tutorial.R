setwd("C:/Abbie/research/seurat/tutorial/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


#### 1. Setup the Seurat Object ####

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(
  counts = pbmc.data, 
  project = "pbmc3k", 
  min.cells = 3, 
  min.features = 200
)
pbmc
## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features, 0 variable features)
##  1 layer present: counts


#### 2. Standard pre-processing workflow ####

### 2a. QC and selecting cells for further analysis ###

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(
  pbmc, 
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
  ncol = 3
)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(
  pbmc, 
  feature1 = "nCount_RNA", 
  feature2 = "percent.mt"
)
plot2 <- FeatureScatter(
  pbmc, 
  feature1 = "nCount_RNA", 
  feature2 = "nFeature_RNA"
)
plot1 + plot2


#### 3. Normalizing the data ####

pbmc <- NormalizeData(
  pbmc, normalization.method = "LogNormalize", 
  scale.factor = 10000
)


#### 4. Identification of highly variable features (feature selection) ####

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#### 5. Scaling the data ####

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#### 6. Perform linear dimensional reduction ####

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca") + NoLegend()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


#### 7. Determine the ‘dimensionality’ of the dataset ####

ElbowPlot(pbmc)


#### 8. Cluster the cells ####

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2638
## Number of edges: 95965
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8723
## Number of communities: 9
## Elapsed time: 0 seconds

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
## AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 
##                2                3                2                1 
## AAACCGTGTATGCG-1 
##                6 
## Levels: 0 1 2 3 4 5 6 7 8


#### 9. Run non-linear dimensional reduction (UMAP/tSNE) ####

pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "output/pbmc_tutorial.rds")


#### 10. Finding differentially expressed features (cluster biomarkers) ####

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)
##             p_val avg_log2FC pct.1 pct.2    p_val_adj
## IL32 2.593535e-91  1.3221171 0.949 0.466 3.556774e-87
## LTB  7.994465e-87  1.3450377 0.981 0.644 1.096361e-82
## CD3D 3.922451e-70  1.0562099 0.922 0.433 5.379250e-66
## IL7R 1.130870e-66  1.4256944 0.748 0.327 1.550876e-62
## LDHB 4.082189e-65  0.9765875 0.953 0.614 5.598314e-61

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)
##                       p_val avg_log2FC pct.1 pct.2     p_val_adj
## FCGR3A        2.150929e-209   6.832372 0.975 0.039 2.949784e-205
## IFITM3        6.103366e-199   6.181000 0.975 0.048 8.370156e-195
## CFD           8.891428e-198   6.052575 0.938 0.037 1.219370e-193
## CD68          2.374425e-194   5.493138 0.926 0.035 3.256286e-190
## RP11-290F20.3 9.308287e-191   6.335402 0.840 0.016 1.276538e-186

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
## # A tibble: 7,046 × 7
## # Groups:   cluster [9]
##        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
##        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
##  1 1.74e-109       1.19 0.897 0.593 2.39e-105 0       LDHB     
##  2 1.17e- 83       2.37 0.435 0.108 1.60e- 79 0       CCR7     
##  3 8.94e- 79       1.09 0.838 0.403 1.23e- 74 0       CD3D     
##  4 3.05e- 53       1.02 0.722 0.399 4.19e- 49 0       CD3E     
##  5 3.28e- 49       2.10 0.333 0.103 4.50e- 45 0       LEF1     
##  6 6.66e- 49       1.25 0.623 0.358 9.13e- 45 0       NOSIP    
##  7 9.31e- 44       2.02 0.328 0.11  1.28e- 39 0       PRKCQ-AS1
##  8 4.69e- 43       1.53 0.435 0.184 6.43e- 39 0       PIK3IP1  
##  9 1.47e- 39       2.70 0.195 0.04  2.01e- 35 0       FHIT     
## 10 2.44e- 33       1.94 0.262 0.087 3.34e- 29 0       MAL      
## # ℹ 7,036 more rows

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


#### 11. Assigning cell type identity to clusters ####

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

saveRDS(pbmc, file = "output/pbmc3k_final.rds")