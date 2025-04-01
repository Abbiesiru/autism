setwd("C:/Abbie/research/seurat/GSE76381/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

seurat_obj <- readRDS("example/GSE76381.mouse_neurons.CPM.rds") # saved seurat object