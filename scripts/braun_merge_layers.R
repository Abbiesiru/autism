base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/braun",
                   "gauss" = "/home/abbiew/single_cell/braun",
                   "."
)

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_combined.rds")
seurat_obj <- readRDS(seurat_obj_path)

library(Seurat)
library(Matrix)

# Get all layer names
layer_names <- Layers(seurat_obj[["RNA"]])

# Filter for counts and data layers
counts_layers <- grep("^counts\\.", layer_names, value = TRUE)
data_layers <- grep("^data\\.", layer_names, value = TRUE)

# Initialize merged matrices
merged_counts <- NULL
merged_data <- NULL

# Merge counts layers
for (layer in counts_layers) {
  cat("Processing", layer, "\n")
  mat <- Seurat:::GetMultiLayerData(seurat_obj, assay = "RNA", slot = "counts", layer = layer)
  merged_counts <- if (is.null(merged_counts)) mat else Matrix::cbind2(merged_counts, mat)
}

# Merge data layers
for (layer in data_layers) {
  cat("Processing", layer, "\n")
  mat <- Seurat:::GetMultiLayerData(seurat_obj, assay = "RNA", slot = "data", layer = layer)
  merged_data <- if (is.null(merged_data)) mat else Matrix::cbind2(merged_data, mat)
}
