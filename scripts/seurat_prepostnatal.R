setwd("C:/Abbie/research/seurat/prepostnatal/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)


#### 1. Setup the Seurat Object ####

# Read matrix
matrix <- readMM("matrix.mtx")

# Read features
features <- read.delim("features.tsv", header = FALSE)
barcodes <- read.delim("barcodes.tsv", header = FALSE)

# Read metadata
metadata <- read.delim("meta.tsv")

# Assign row and column names to matrix
rownames(matrix) <- make.unique(features$V1) 
colnames(matrix) <- barcodes$V1

# Create Seurat obj
seurat_obj <- CreateSeuratObject(
  counts = matrix,
  meta.data = metadata,
  project = "PrePostNatal"
)


#### 3. Normalizing the data ####

seurat_obj <- NormalizeData(
  seurat_obj, normalization.method = "LogNormalize", 
  scale.factor = 10000
)


#### 4. Identification of highly variable features (feature selection) ####

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)


#### 5. Scaling the data ####

all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)


#### 6. Visualize gene expression by developmental stage ####

# Read risk genes
risk_genes <- read.csv("autism_risk_genes_combined.csv")

FeaturePlot(seurat_obj, features = risk_genes, split.by = "Age")
VlnPlot(seurat_obj, features = risk_genes, group.by = "Age", pt.size = 0)

# Read umap coords
tmp <- read.delim("UMAP.coords.tsv.gz", sep = "\t", header = FALSE)
umap <- tmp[,-1]
rownames(umap) <- tmp[,1]
colnames(umap) <- c("UMAP_1", "UMAP_2")

seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = "RNA")
