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
                                          cell_types, cell_col, age_col, dataset_label) {
  # Transpose matrices: now rows = cells, columns = genes
  rank_df <- as.data.frame(t(rank_mat))
  expr_df <- as.data.frame(t(expr_mat))
  
  # Add cell names
  rank_df$Cell <- rownames(rank_df)
  expr_df$Cell <- rownames(expr_df)
  
  # Merge expression and rank by Cell
  merged <- merge(rank_df, expr_df, by = "Cell", suffixes = c("_rank", "_expr"))
  
  # Extract metadata from Seurat
  meta <- seurat_obj@meta.data
  meta$Cell <- rownames(meta)
  
  # Filter metadata to desired cell types
  meta <- meta[meta[[cell_col]] %in% cell_types, c("Cell", cell_col, age_col)]
  colnames(meta) <- c("Cell", "Cell_Type", "Age")
  
  # Merge all together
  merged <- merge(merged, meta, by = "Cell")
  
  # Get list of genes
  gene_names <- sub("_rank$", "", grep("_rank$", colnames(merged), value = TRUE))
  
  # Initialize output storage
  summary_list <- list()
  
  for (gene in gene_names) {
    rank_col <- paste0(gene, "_rank")
    expr_col <- paste0(gene, "_expr")
    
    df <- merged[, c("Cell_Type", "Age", rank_col, expr_col)]
    colnames(df) <- c("Cell_Type", "Age", "Rank", "Expr")
    
    # Aggregate by Cell_Type and Age
    agg_df <- aggregate(df[, c("Rank", "Expr")],
                        by = list(Gene = gene, Cell_Type = df$Cell_Type, Age = df$Age),
                        FUN = function(x) c(mean = mean(x, na.rm = TRUE), pct = mean(x > 0, na.rm = TRUE) * 100))
    
    # Unpack list-columns into separate columns
    agg_df$avg_rank <- sapply(agg_df$Rank, function(x) x[1])
    agg_df$percent_expressing <- sapply(agg_df$Expr, function(x) x[2])
    agg_df$Dataset <- dataset_label
    
    # Drop intermediate list columns
    agg_df$Rank <- NULL
    agg_df$Expr <- NULL
    
    summary_list[[gene]] <- agg_df
  }
  
  summary_all <- do.call(rbind, summary_list)
  rownames(summary_all) <- NULL
  return(summary_all)
}


# Neuronal lineage types
neuronal_types_b <- c("Neuroblast", "Neuron", "Neuronal IPC")
neuronal_types_v <- c("ExNeu", "IN")

# Process both datasets
neuronal_b <- process_rank_and_percent(
  rank_mat = rank_mat_b,
  expr_mat = expr_mat_b,
  seurat_obj = b,
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
