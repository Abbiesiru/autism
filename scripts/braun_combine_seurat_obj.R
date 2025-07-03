setwd("/home/abbiew/single_cell/braun")
seurat_paths <- list.files("./5_Saved_seurat_obj", pattern = "*.rds", full.names = TRUE)
tmp <- read.csv("/home/abbiew/single_cell/autism_risk_genes_combined.csv", sep = ",", header = TRUE)
asd_risk_genes <- tmp$Gene

library(Seurat)

process_and_subset <- function(path, genes) {
  seurat_obj <- readRDS(path)
  
  seurat_obj@assays$RNA@layers$data@x <- log2(seurat_obj@assays$RNA@layers$data@x + 1)
  Idents(seurat_obj, cells = seurat_obj$Cell_name) <- seurat_obj$Cell_type_raw
  
  temp_data <- as.matrix(seurat_obj@meta.data[, c("tSNE_x", "tSNE_y")])
  colnames(temp_data) <- c("tSNE_1", "tSNE_2");
  seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = temp_data, key = "tSNE_",
                                                    assay = DefaultAssay(seurat_obj))

  seurat_obj <- seurat_obj[asd_risk_genes, ]

  return(seurat_obj)
}

processed_list <- lapply(seurat_paths, process_and_subset, genes = asd_risk_genes)
combined_seurat <- merge(x = processed_list[[1]], y = processed_list[-1])

saveRDS(combined_seurat, "seurat_obj_subset_combined.rds")
