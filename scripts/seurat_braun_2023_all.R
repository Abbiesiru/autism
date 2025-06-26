setwd("C:/Abbie/research/seurat/braun_2023/cellranger_output_split")

library(Seurat)
library(Matrix)
library(ggplot2)

data_dirs <- list(
  "brain" = "Brain/matrix.mtx",
  "cerebellum" = "Cerebellum/matrix.mtx",
  "diencephalon" = "Diencephalon/matrix.mtx",
  "forebrain" = "Forebrain/matrix.mtx",
  "Head"= "Head/matrix.mtx",
  "Hindbrain" = "Hindbrain/matrix.mtx"
)
