base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat",
                   "gauss" = "/home/abbiew/single_cell",
                   "."
)

output_folder <- "plots"
output_dir <- file.path(base_dir, output_folder)

# Create folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)

seurat_obj_path <- file.path(base_dir, "braun/seurat_obj_merged_layers_only_processed.rds")
b <- readRDS(seurat_obj_path)
seurat_obj_path <- file.path(base_dir, "velmeshev/seurat_obj_subset_common_genes.rds")
v <- readRDS(seurat_obj_path)
rank_mat_b <- readRDS(file.path(base_dir, "/braun/cell_rankings_braun.rds"))
rank_mat_v <- readRDS(file.path(base_dir, "/prepostnatal/cell_rankings_velmeshev.rds"))
expr_mat_b <- GetAssayData(b, assay = "merged", slot = "data")
expr_mat_v <- GetAssayData(v, assay = "RNA", slot = "data")


# Helper function to process both rank and percent expressed
process_rank_and_percent <- function(rank_mat, expr_mat, seurat_obj,
                                     gene_col, age_col, cell_col, cell_types, dataset_label) {
  rank_df <- as.data.frame(t(rank_mat))  # cells × genes
  rank_df$Cell <- rownames(rank_df)
  
  expr_df <- as.data.frame(t(expr_mat))  # cells × genes
  expr_df$Cell <- rownames(expr_df)
  
  # Join rank + expression
  merged <- left_join(rank_df, expr_df, by = "Cell", suffix = c("_rank", "_expr"))
  
  # Add metadata
  meta <- seurat_obj@meta.data %>%
    rownames_to_column(var = "Cell") %>%  # cell barcodes become a column named "Cell"
    select(Cell, !!sym(cell_col), !!sym(age_col)) %>%
    filter(!!sym(cell_col) %in% cell_types)

  
  merged <- left_join(merged, meta, by = "Cell") %>%
    drop_na()
  
  # Pivot rank and expr side-by-side
  rank_long <- merged %>%
    pivot_longer(cols = ends_with("_rank"), names_to = "Gene_rank", values_to = "Rank") %>%
    mutate(Gene = sub("_rank$", "", Gene_rank)) %>%
    select(-Gene_rank)
  
  expr_long <- merged %>%
    pivot_longer(cols = ends_with("_expr"), names_to = "Gene_expr", values_to = "Expr") %>%
    mutate(Gene = sub("_expr$", "", Gene_expr)) %>%
    select(-Gene_expr)
  
  # Merge rank + expr
  long_df <- left_join(rank_long, expr_long, by = c("Cell", "Gene", "Cell", "Age_Range" = age_col, "Cell_type" = cell_col))
  
  # Compute summaries
  summary <- long_df %>%
    group_by(Gene, !!sym(cell_col), !!sym(age_col)) %>%
    summarize(
      avg_rank = mean(Rank, na.rm = TRUE),
      percent_expressing = mean(Expr > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    rename(Cell_Type = !!sym(cell_col), Age = !!sym(age_col)) %>%
    mutate(Dataset = dataset_label)
  
  return(summary)
}

# Neuronal lineage types
neuronal_types_b <- c("Neuroblast", "Neuron", "Neuronal IPC")
neuronal_types_v <- c("ExNeu", "IN")

# Process both datasets
neuronal_b <- process_rank_and_percent(
  rank_mat = rank_mat_b,
  expr_mat = expr_mat_b,
  seurat_obj = b,
  gene_col = "Gene",
  age_col = "Developmental_week",
  cell_col = "Cell_type_raw",
  cell_types = neuronal_types_b,
  dataset_label = "Braun"
)

neuronal_v <- process_rank_and_percent(
  rank_mat = rank_mat_v,
  expr_mat = expr_mat_v,
  seurat_obj = v,
  gene_col = "Gene",
  age_col = "Age_Range",
  cell_col = "Lineage",
  cell_types = neuronal_types_v,
  dataset_label = "Velmeshev"
)

# Combine and set factor levels for age
age_levels_combined <- c(
  "5", "5.5", "6", "6.6", "6.7", "6.9", "7", "7.5", "8", "8.1", "8.5", "9.2", "9.5", "10", "11.5", "12", "13", "14",
  "2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult"
)

combined_rank <- bind_rows(neuronal_b, neuronal_v) %>%
  mutate(Age = factor(Age, levels = age_levels_combined))

# Final plot
ggplot(combined_rank, aes(x = Age, y = avg_rank, color = Dataset, size = percent_expressing)) +
  geom_point(alpha = 0.8) +
  geom_line(aes(group = interaction(Dataset, Cell_Type)), linewidth = 1) +
  facet_wrap(~Gene + Cell_Type, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_size_continuous(name = "% Expressing", range = c(1, 6)) +
  labs(
    title = "Expression Rank Over Development (Neuronal Lineages)",
    x = "Developmental Age",
    y = "Scaled Rank (0–1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )
