library(Seurat)
library(Matrix)

data_dir <- "C:/Abbie/research/seurat/braun_2023/filtered_10x_output_5"

filtered_data <- ReadMtx(
  mtx = file.path(data_dir, "matrix.mtx"),
  features = file.path(data_dir, "features.tsv"),
  cells = file.path(data_dir, "barcodes.tsv"),
  feature.column = 2
)

seurat_obj <- CreateSeuratObject(counts = filtered_data)

# Get raw counts matrix (sparse matrix)
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# Get normalized data
norm_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

# Get scaled data (used for PCA, clustering, etc.)
scaled_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "scale.data")
