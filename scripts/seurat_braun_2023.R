setwd("C:/Abbie/research/seurat/braun_2023/")

library(Seurat)
library(Matrix)
library(ggplot2)

# Load the sparse matrix, genes & barcodes
expr_matrix <- readMM("cellranger_output_2/matrix.mtx") 
expr_matrix <- t(expr_matrix)
genes <- read.delim("cellranger_output_2/features.tsv", header = FALSE)
barcodes <- read.delim("cellranger_output_2/barcodes.tsv", header = FALSE)

rownames(expr_matrix) <- genes$V2   # set gene names as row names
colnames(expr_matrix) <- barcodes$V1 # set cell barcodes as col names

seurat_obj <- CreateSeuratObject(counts = expr_matrix)

cell_anno <- read.delim("cell_anno.txt", header = TRUE, sep = "\t")

# Set rownames as barcodes to align with Seurat cells
rownames(cell_anno) <- cell_anno$Cell_name

# Subset to cells in Seurat object
cell_anno <- cell_anno[colnames(seurat_obj), , drop = FALSE]

# Add annotation as metadata
seurat_obj <- AddMetaData(seurat_obj, metadata = cell_anno)

# Extract the t-SNE coords from metadata as a matrix
colnames(cell_anno)[colnames(cell_anno) == "tSNE_x"] <- "tSNE_1"
colnames(cell_anno)[colnames(cell_anno) == "tSNE_y"] <- "tSNE_2"
tsne_coords <- as.matrix(cell_anno[, c("tSNE_1", "tSNE_2")])

# Create a new DimReduc object
seurat_obj[["tSNE_plot"]] <- CreateDimReducObject(
  embeddings = tsne_coords,
  key = "TSNE_",
  assay = DefaultAssay(seurat_obj)
)

ggplot(cell_anno, aes(x = tSNE_1, y = tSNE_2, color = Region)) +
  geom_point() +
  theme_minimal()

