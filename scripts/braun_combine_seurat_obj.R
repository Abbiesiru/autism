setwd("/home/abbiew/single_cell/braun")
seurat_paths <- list.files("./5_Saved_seurat_obj", pattern = "*.rds", full.names = TRUE)
tmp <- read.csv("/home/abbiew/single_cell/autism_risk_genes_combined.csv", sep = ",", header = TRUE)
asd_risk_genes <- tmp$Gene

library(Seurat)
library(AUCell)

cell_rankings_list <- list()

process_and_subset <- function(path, genes) {
  seurat_obj <- readRDS(path)

  # 1. rank the genes in each cell, subset & export
  expr_mat <- GetAssayData(object = seurat_obj, assay = "RNA", slot = "counts")
  cell_rankings <- AUCell_buildRankings(
    expr_mat,
    plotStats = FALSE,
    keepZeroesAsNA = FALSE,
    verbose = TRUE
  )
  subset_rankings <- cell_rankings[genes, , drop = FALSE]
  cell_rankings_list[[basename(path)]] <<- subset_rankings 
  

  # 2. log2 transformation
  seurat_obj@assays$RNA@layers$data@x <- log2(seurat_obj@assays$RNA@layers$data@x + 1)
  Idents(seurat_obj, cells = seurat_obj$Cell_name) <- seurat_obj$Cell_type_raw

  # 3. add tsne coordinates
  temp_data <- as.matrix(seurat_obj@meta.data[, c("tSNE_x", "tSNE_y")])
  colnames(temp_data) <- c("tSNE_1", "tSNE_2");
  seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = temp_data, key = "tSNE_",
                                                    assay = DefaultAssay(seurat_obj))

  # 4. subset to ASD risk genes only 
  seurat_obj <- seurat_obj[asd_risk_genes, ]

  return(seurat_obj)
}

processed_list <- lapply(seurat_paths, process_and_subset, genes = asd_risk_genes)
combined_seurat <- merge(x = processed_list[[1]], y = processed_list[-1])

merged_cell_rankings <- do.call(cbind, cell_rankings_list)
saveRDS(merged_cell_rankings, file = "cell_rankings_combined.rds")

saveRDS(combined_seurat, "seurat_obj_subset_combined.rds")
