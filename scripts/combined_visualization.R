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
library(Matrix)

genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")

seurat_obj_path <- file.path(base_dir, "braun/seurat_obj_merged_layers_only_processed.rds")
b <- readRDS(seurat_obj_path)
seurat_obj_path <- file.path(base_dir, "velmeshev/seurat_obj_subset_common_genes.rds")
v <- readRDS(seurat_obj_path)
rank_mat_b <- readRDS(file.path(base_dir, "/braun/cell_rankings_braun.rds"))
rank_mat_v <- readRDS(file.path(base_dir, "/velmeshev/cell_rankings_velmeshev.rds"))
expr_mat_b <- GetAssayData(b, assay = "merged", slot = "data")
expr_mat_v <- GetAssayData(v, assay = "RNA", slot = "data")


# Helper function to process both rank and percent expressed
process_rank_and_percent <- function(rank_mat, expr_mat, seurat_obj,
                                             cell_types, cell_col, age_col, dataset_label) {
  # Step 1: Get and filter metadata
  meta <- seurat_obj@meta.data
  meta$Cell <- rownames(meta)
  meta <- meta[meta[[cell_col]] %in% cell_types, c("Cell", cell_col, age_col)]
  colnames(meta) <- c("Cell", "Cell_Type", "Age")
  cells_keep <- meta$Cell
  
  # Step 2: Subset matrices to relevant cells
  rank_mat_sub <- rank_mat[, intersect(colnames(rank_mat), cells_keep), drop = FALSE]
  expr_mat_sub <- expr_mat[, intersect(colnames(expr_mat), cells_keep), drop = FALSE]

  # Step 3: Find common genes
  gene_names <- intersect(rownames(rank_mat_sub), rownames(expr_mat_sub))
  
  # Step 4: Prepare output
  summary_list <- list()
  
  for (gene in gene_names) {
    # Extract rank and expr vectors (1 gene across all cells)
    rank_vec <- rank_mat_sub[gene, ]
    expr_vec <- expr_mat_sub[gene, ]
    
    # Force named numeric vectors
    rank_vec <- as.numeric(rank_vec)
    expr_vec <- as.numeric(expr_vec)
    cells <- colnames(rank_mat_sub)
    
    # Step 5: Build dataframe
    df <- data.frame(
      Cell = cells,
      Rank = rank_vec,
      Expr = expr_vec,
      stringsAsFactors = FALSE
    )
    
    # Join with metadata
    df <- merge(df, meta, by = "Cell")
    
    # Drop missing
    df <- df[complete.cases(df), ]
    if (nrow(df) == 0) next
    
    # Step 6: Check lengths before aggregation
    stopifnot(length(df$Rank) == length(df$Cell_Type))
    stopifnot(length(df$Expr) == length(df$Age))
    
    # Step 7: Aggregate safely
    agg_rank <- aggregate(Rank ~ Cell_Type + Age, data = df, FUN = mean)
    agg_expr <- aggregate(Expr ~ Cell_Type + Age, data = df, FUN = function(x) mean(x > 0) * 100)
    
    # Merge summaries
    agg_df <- merge(agg_rank, agg_expr, by = c("Cell_Type", "Age"))
    agg_df$Gene <- gene
    agg_df$Dataset <- dataset_label
    colnames(agg_df) <- c("Cell_Type", "Age", "avg_rank", "percent_expressing", "Gene", "Dataset")
    
    summary_list[[gene]] <- agg_df
  }
  
  # Combine everything
  if (length(summary_list) == 0) {
    warning("No gene summaries generated.")
    return(data.frame())
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
  age_col = "Age_Range",
  cell_col = "Lineage",
  cell_types = neuronal_types_v,
  dataset_label = "Velmeshev"
)

# Combine both sets of neuronal types
neuronal_data <- combined_rank[combined_rank$Cell_Type %in% c(neuronal_types_b, neuronal_types_v), ]

# Ensure Age is an ordered factor
age_levels <- c(
  "5", "5.5", "6", "6.6", "6.7", "6.9", "7", "7.5", "8", "8.1", 
  "8.5", "9.2", "9.5", "10", "11.5", "12", "13", "14",
  "2nd trimester", "3rd trimester", "0-1 years", "1-2 years",
  "2-4 years", "4-10 years", "10-20 years", "Adult"
)
neuronal_data$Age <- factor(neuronal_data$Age, levels = age_levels, ordered = TRUE)

# --- Loop for each gene ---
for (g in genes_of_interest) {
  # Subset for this gene
  gene_data <- neuronal_data[neuronal_data$Gene == g, ]

  # Aggregate by Gene + Age + Dataset
  agg_rank <- aggregate(
    avg_rank ~ Gene + Age + Dataset,
    data = gene_data,
    FUN = mean,
    na.rm = TRUE
  )
  
  agg_pct <- aggregate(
    percent_expressing ~ Gene + Age + Dataset,
    data = gene_data,
    FUN = mean,
    na.rm = TRUE
  )
  
  # Merge results
  agg_df <- merge(agg_rank, agg_pct, by = c("Gene", "Age", "Dataset"))

  # Plot
  p <- ggplot(agg_df, aes(x = Age, y = avg_rank, color = Dataset, size = percent_expressing)) +
    geom_point(alpha = 0.9) +
    geom_line(aes(group = Dataset), linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_size_continuous(name = "% Expressing", range = c(1, 6)) +
    labs(
      title = paste("Avg Neuronal Expression Rank -", g),
      x = "Developmental Age",
      y = "Scaled Avg Rank (0–1)",
      color = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  # Save
  ggsave(
    filename = file.path(output_dir, paste0("rank_plot_neuronal_", g, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
}
