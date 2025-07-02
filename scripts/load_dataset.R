filename = "/home/abbiew/single_cell/5_Saved_seurat_obj/Braun_Science_2023.scRNA.Head.CPM.rds"
library(Seurat)
seurat_obj <- readRDS(filename)
seurat_obj@assays$RNA@layers$data@x <- log2(seurat_obj@assays$RNA@layers$data@x + 1)
Idents(seurat_obj, cells = seurat_obj$Cell_name) <- seurat_obj$Cell_type_raw;

temp_data <- as.matrix(seurat_obj@meta.data[, c("tSNE_x", "tSNE_y")]);
colnames(temp_data) <- c("tSNE_1", "tSNE_2");
seurat_obj[["tsne"]] <- CreateDimReducObject(embeddings = temp_data, key = "tSNE_",
                                                    assay = DefaultAssay(seurat_obj));

FeaturePlot(seurat_obj, "SORCS1")