library(Seurat)
library(Matrix)

base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/braun",
                   "gauss" = "/home/abbiew/single_cell/braun",
                   ".")

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_combined.rds")
seurat_obj <- readRDS(seurat_obj_path)

rna_assay <- seurat_obj[["RNA"]]

layer_names <- Layers(rna_assay)

counts_layers <- grep("^counts\\.", layer_names, value = TRUE)
data_layers <- grep("^data\\.", layer_names, value = TRUE)

# Merge counts matrices
merged_counts <- NULL
for (layer in counts_layers) {
  cat("Processing counts layer:", layer, "\n")
  mat <- rna_assay@layers[[layer]]
  
  # Ensure colnames exist and are unique with prefix
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0(layer, "_cell", seq_len(ncol(mat)))
  } else {
    colnames(mat) <- paste(layer, colnames(mat), sep = "_")
  }
  
  # Ensure rownames exist (use RNA assay gene names as fallback)
  if (is.null(rownames(mat))) {
    rownames(mat) <- rownames(rna_assay)
  }
  
  merged_counts <- if (is.null(merged_counts)) mat else Matrix::cbind2(merged_counts, mat)
}

# Merge data matrices
merged_data <- NULL
for (layer in data_layers) {
  cat("Processing data layer:", layer, "\n")
  mat <- rna_assay@layers[[layer]]
  
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0(layer, "_cell", seq_len(ncol(mat)))
  } else {
    colnames(mat) <- paste(layer, colnames(mat), sep = "_")
  }
  
  if (is.null(rownames(mat))) {
    rownames(mat) <- rownames(rna_assay)
  }
  
  merged_data <- if (is.null(merged_data)) mat else Matrix::cbind2(merged_data, mat)
}

original_cells <- colnames(seurat_obj[["RNA"]])
colnames(merged_counts) <- original_cells
colnames(merged_data) <- original_cells

# Create new assay with merged counts
merged_assay <- CreateAssayObject(counts = merged_counts)


if (!is.null(merged_data)) {
  merged_assay@data <- as(merged_data, "CsparseMatrix")
}

# Add new assay to Seurat object
seurat_obj[["merged"]] <- merged_assay
DefaultAssay(seurat_obj) <- "merged"
seurat_obj[["RNA"]] <- NULL

# Save updated Seurat object
saveRDS(seurat_obj, file = file.path(base_dir, "seurat_obj_merged_layers.rds"))


