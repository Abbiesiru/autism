base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/prepostnatal",
                   "gauss" = "/home/abbiew/single_cell/velmeshev",
                   "."
)

output_folder <- "plots"
output_dir <- file.path(base_dir, output_folder)

# Create folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(Seurat)
library(ggplot2)

setwd("C:/Abbie/research/seurat")
asd_risk_genes <- read.csv("./autism_risk_genes_combined.csv", sep = ",", header = TRUE)

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_common_genes.rds")
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj$Lineage <- recode(
  seurat_obj$Lineage,
  "AST" = "Astrocytes",
  "ExNeu" = "Excitatory neurons",
  "GLIALPROG" = "Glial progenitors",
  "IN" = "Inhibitory neurons",
  "MG" = "Microglia",
  "OL" = "Oligodendrocytes",
  "OPC" = "Oligodendrocyte precursors",
  "OUT" = "Outliers",
  "VASC" = "Vascular cells"
)

subset_obj <- subset(seurat_obj, subset = Age_Range == "1-2 years")
table(subset_obj$Lineage)


df <- as.data.frame(table(subset_obj$Lineage))
colnames(df) <- c("CellType", "Count")

ggplot(df, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Cell Type Distribution (Age 1–2 Years)",
       x = "Cell Type", y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


df <- as.data.frame(table(seurat_obj$Lineage))
colnames(df) <- c("CellType", "Count")

ggplot(df, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Cell Type Distribution (Age 1–2 Years)",
       x = "Cell Type", y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DimPlot(
  subset_obj,
  reduction = "umap",
  group.by = "Region_Broad",     # color by cell type
  label = TRUE,               # add labels to clusters
  label.size = 3,
  pt.size = 1.5
  )

subset_astro_age <- subset(
  seurat_obj,
  subset = Age_Range == "1-2 years" & Lineage == "Astrocytes"
)

subset_astro <- subset(
  seurat_obj,
  Lineage == "Astrocytes"
)
  
  
 p <- DotPlot(
    object = subset_astro,
    features = "SORCS2",
    group.by = "Region_Broad",
    cols = c("grey", "blue")
  ) +
    labs(x = "Gene", y = "Region") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("Regional Expresson of SORCS2 in Astrocytes")
  

  ggsave(
    filename = file.path(output_dir, "dot_plot_SORCS2_astro_region_broad.pdf"),
    plot = p,
    width = 10,
    height = 10
  )

  df <- as.data.frame(table(subset_astro$Region_Broad))
  colnames(df) <- c("Region", "Count")
  
  p <- ggplot(df, aes(x = reorder(Region, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "forestgreen") +
    theme_minimal() +
    labs(title = "Astrocyte Count by Brain Region",
         x = "Brain Region", y = "Number of Astrocytes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(output_dir, "bar_chart_astro_region_broad.pdf"),
    plot = p,
    width = 10,
    height = 10
  )
  
  
  library(dplyr)
  
  
  age_order <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")
  
  # Set Age_Range as factor with correct order
  seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_order)
  
  
  donor_counts <- seurat_obj@meta.data %>%
    distinct(Individual, Age_Range) %>%
    count(Age_Range, name = "Num_Donors")
  
  p <- ggplot(donor_counts, aes(x = Age_Range, y = Num_Donors)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_minimal() +
    labs(title = "Number of Donors per Age Range",
         x = "Age Range", y = "Number of Donors") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(output_dir, "donors_per_age_range.pdf"),
    plot = p,
    width = 10,
    height = 5
  )
  
  sample_counts <- seurat_obj@meta.data %>%
    distinct(Sample, Age_Range) %>%
    count(Age_Range, name = "Num_Samples")
  
  p <- ggplot(sample_counts, aes(x = Age_Range, y = Num_Samples)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    theme_minimal() +
    labs(title = "Number of Samples per Age Range",
         x = "Age Range", y = "Number of Samples") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(output_dir, "samples_per_age_range.pdf"),
    plot = p,
    width = 10,
    height = 5
  )