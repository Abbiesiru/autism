setwd("C:/Abbie/research/seurat/braun_2023/")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SeuratDisk)
library(zellkonverter)



#### 1. Setup the Seurat Object ####

# Load file
Convert("human_dev_GRCh38-3.0.0.h5ad", dest = "h5seurat", overwrite = TRUE)
zellkonverter::readH5AD("human_dev.h5ad", use_hdf5 = FALSE)
seurat_obj <- Read10X_h5("human_dev_GRCh38-3.0.0.h5seurat")
