setwd("/home/abbiew/single_cell/braun")
seurat_paths <- list.files("./5_Saved_seurat_obj", pattern = "*.rds", full.names = TRUE)
tmp <- read.csv("/home/abbiew/single_cell/autism_risk_genes_combined.csv", sep = ",", header = TRUE)
asd_risk_genes <- tmp$Gene

library(Seurat)

process_and_subset <- function(path, genes) {
  seurat_obj <- readRDS(path)

  # 1. log2 transformation
  expr_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  log_expr <- log2(expr_data + 1)
  seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", slot = "data", new.data = log_expr)

  Idents(seurat_obj, cells = seurat_obj$Cell_name) <- seurat_obj$Cell_type_raw

  # 2. add tsne coordinates
  temp_data <- as.matrix(seurat_obj@meta.data[, c("tSNE_x", "tSNE_y")])
  colnames(temp_data) <- c("tSNE_1", "tSNE_2");
  seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = temp_data, key = "tSNE_",
                                                    assay = DefaultAssay(seurat_obj))

  # 3. subset to ASD risk genes only 
  seurat_obj <- seurat_obj[rownames(seurat_obj) %in% asd_risk_genes, ]

  return(seurat_obj)
}

processed_list <- list()
for (i in seq_along(seurat_paths)) {
  processed_list[[i]] <- process_and_subset(seurat_paths[[i]], genes = asd_risk_genes)
  gc()  # Trigger garbage collection
}

combined_seurat <- Reduce(function(x, y) merge(x, y), processed_list)

saveRDS(combined_seurat, "seurat_obj_subset_combined.rds")
