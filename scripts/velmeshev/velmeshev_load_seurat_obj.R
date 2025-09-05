base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/prepostnatal",
                   "gauss" = "/home/abbiew/single_cell/velmeshev",
                   "."
)
setwd("C:/Abbie/research/seurat/prepostnatal")

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

mat <- readMM("matrix.mtx.gz")
barcodes <- read.delim("barcodes.tsv.gz", header = FALSE)
features <- read.delim("features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(mat) <- make.unique(features$V1) # gene names
colnames(mat) <- barcodes$V1 # cell barcodes

seurat_obj <- CreateSeuratObject(counts = mat, project = "prepostnatal", meta.data=meta)

umap <- read.delim("UMAP.coords.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
colnames(umap) <- c("barcode", "UMAP_1", "UMAP_2")
rownames(umap) <- umap$barcode
umap$barcode <- NULL
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = "RNA")

seurat_obj <- NormalizeData(seurat_obj)

saveRDS(seurat_obj, "seurat_obj.rds")

### subset seurat_obj ###

ASD_risk_genes <- read.csv("autism_risk_genes_combined.csv", sep = ",", header = TRUE)

seurat_obj_path <- file.path(base_dir, "seurat_obj.rds")
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj_subset <- seurat_obj[ASD_risk_genes$Gene, ]

# sort age ranges from earliest to latest
age_levels <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")
seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_levels)

# saveRDS(seurat_obj_subset, "seurat_obj_subset.rds")
# seurat_obj_path <- file.path(base_dir, "seurat_obj_subset.rds")
# seurat_obj_subset <- readRDS(seurat_obj_path)

# subset to 266 common asd risk genes
common_asd_risk_genes <- readRDS("C:/Abbie/research/seurat/common_asd_risk_genes.rds")
seurat_obj_subset <- seurat_obj_subset[common_asd_risk_genes, ]
saveRDS(seurat_obj_subset, file.path(base_dir, "seurat_obj_subset_common_genes.rds"))
