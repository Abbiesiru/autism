setwd("/Users/abbiesiru/Desktop/research/seurat/GSE104276/")

library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)
library(Matrix)
library(ggplot2)

example <- readRDS("GSE104276.CPM_example.rds") # saved seurat object
tmp1 <- read_excel("GSE104276_all_pfc_2394_UMI_count_NOERCC.xlsx")
tmp2 <- read.delim("cell_info.GSE104276.txt", header = TRUE, sep = "\t")
matrix <- tmp1[, colnames(tmp1) %in% rownames(tmp2)]
metadata <- tmp2[colnames(matrix), ]

seurat_obj <- CreateSeuratObject(
  counts = as.matrix(matrix),
  meta.data = metadata,
  project = "GSE104276"
)


