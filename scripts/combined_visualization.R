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

age_levels <- c(
  "5", "5.5", "6", "6.6", "6.7", "6.9", "7", "7.5", "8", "8.1", 
  "8.5", "9.2", "9.5", "10", "11.5", "12", "13", "14",
  "2nd trimester", "3rd trimester", "0-1 years", "1-2 years",
  "2-4 years", "4-10 years", "10-20 years", "Adult"
)

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


#### 1. Neurons ####

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
     scale_size_continuous(name = "% Expressing", limits = c(0, 100), range = c(1, 6)) +
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


#### 2. OPCs ####

# OPC for Braun: Subclass == "OPC"
opc_b <- process_rank_and_percent(
  rank_mat = rank_mat_b,
  expr_mat = expr_mat_b,
  seurat_obj = b,
  age_col = "Developmental_week",
  cell_col = "Subclass",
  cell_types = "OPC",
  dataset_label = "Braun"
)

# OPC for Velmeshev: Lineage == "OPC"
opc_v <- process_rank_and_percent(
  rank_mat = rank_mat_v,
  expr_mat = expr_mat_v,
  seurat_obj = v,
  age_col = "Age_Range",
  cell_col = "Lineage",
  cell_types = "OPC",
  dataset_label = "Velmeshev"
)

# Combine Braun + Velmeshev OPC summaries
opc_data <- rbind(opc_b, opc_v)

# Ensure Age is ordered factor
opc_data$Age <- factor(opc_data$Age, levels = age_levels, ordered = TRUE)

# --- Plot loop (same as before) ---
for (g in genes_of_interest) {
  gene_data <- opc_data[opc_data$Gene == g, ]
  
  # Aggregate by Gene + Age + Dataset (technically redundant since only 1 cell type per dataset)
  agg_rank <- aggregate(avg_rank ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  agg_pct  <- aggregate(percent_expressing ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  
  agg_df <- merge(agg_rank, agg_pct, by = c("Gene", "Age", "Dataset"))
  
  # Plot
  p <- ggplot(agg_df, aes(x = Age, y = avg_rank, color = Dataset, size = percent_expressing)) +
    geom_point(alpha = 0.9) +
    geom_line(aes(group = Dataset), linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_size_continuous(name = "% Expressing", limits = c(0, 100), range = c(1, 6)) +
    labs(
      title = paste("Avg OPC Expression Rank -", g),
      x = "Developmental Age",
      y = "Scaled Avg Rank (0–1)",
      color = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  # Save to output_dir
  ggsave(
    filename = file.path(output_dir, paste0("rank_plot_OPC_", g, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
}


#### 3. non-neuronal ####

# Get all cell types present in datasets
all_b_types <- unique(b@meta.data$Cell_type_raw)
all_v_types <- unique(v@meta.data$Lineage)

# Define non-neuronal types (all except neuronal)
non_neuronal_b <- setdiff(all_b_types, neuronal_types_b)
non_neuronal_v <- setdiff(all_v_types, neuronal_types_v)

nonneuronal_b <- process_rank_and_percent(
  rank_mat = rank_mat_b,
  expr_mat = expr_mat_b,
  seurat_obj = b,
  age_col = "Developmental_week",
  cell_col = "Cell_type_raw",
  cell_types = non_neuronal_b,
  dataset_label = "Braun"
)

# Process Velmeshev non-neuronal
nonneuronal_v <- process_rank_and_percent(
  rank_mat = rank_mat_v,
  expr_mat = expr_mat_v,
  seurat_obj = v,
  age_col = "Age_Range",
  cell_col = "Lineage",
  cell_types = non_neuronal_v,
  dataset_label = "Velmeshev"
)

# Combine datasets
nonneuronal_data <- rbind(nonneuronal_b, nonneuronal_v)

# Ensure Age is an ordered factor (reuse age_levels vector)
nonneuronal_data$Age <- factor(nonneuronal_data$Age, levels = age_levels, ordered = TRUE)

# --- Plotting loop for non-neuronal ---
for (g in genes_of_interest) {
  gene_data <- nonneuronal_data[nonneuronal_data$Gene == g, ]
  
  # Aggregate by Gene + Age + Dataset
  agg_rank <- aggregate(avg_rank ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  agg_pct <- aggregate(percent_expressing ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  
  agg_df <- merge(agg_rank, agg_pct, by = c("Gene", "Age", "Dataset"))
  
  # Plot
  p <- ggplot(agg_df, aes(x = Age, y = avg_rank, color = Dataset, size = percent_expressing)) +
    geom_point(alpha = 0.9) +
    geom_line(aes(group = Dataset), linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_size_continuous(name = "% Expressing", limits = c(0, 100), range = c(1, 6)) +
    labs(
      title = paste("Avg Non-Neuronal Expression Rank -", g),
      x = "Developmental Age",
      y = "Scaled Avg Rank (0–1)",
      color = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0("rank_plot_nonneuronal_", g, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
}

#### 4. Pre-astrocytes & astrocytes ####

# Pre-astrocytes for Braun: Subclass == "PREAC"
ast_b <- process_rank_and_percent(
  rank_mat = rank_mat_b,
  expr_mat = expr_mat_b,
  seurat_obj = b,
  age_col = "Developmental_week",
  cell_col = "Subclass",
  cell_types = "PREAC",
  dataset_label = "Braun"
)

# Astrocytes for Velmeshev: Lineage == "AST"
ast_v <- process_rank_and_percent(
  rank_mat = rank_mat_v,
  expr_mat = expr_mat_v,
  seurat_obj = v,
  age_col = "Age_Range",
  cell_col = "Lineage",
  cell_types = "AST",
  dataset_label = "Velmeshev"
)

# Combine Braun + Velmeshev astrocyte summaries
ast_data <- rbind(ast_b, ast_v)

# Ensure Age is ordered factor
ast_data$Age <- factor(ast_data$Age, levels = age_levels, ordered = TRUE)

# --- Plot loop (same as before) ---
for (g in genes_of_interest) {
  gene_data <- ast_data[ast_data$Gene == g, ]
  
  # Aggregate by Gene + Age + Dataset (technically redundant since only 1 cell type per dataset)
  agg_rank <- aggregate(avg_rank ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  agg_pct  <- aggregate(percent_expressing ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  
  agg_df <- merge(agg_rank, agg_pct, by = c("Gene", "Age", "Dataset"))
  
  # Plot
  p <- ggplot(agg_df, aes(x = Age, y = avg_rank, color = Dataset, size = percent_expressing)) +
    geom_point(alpha = 0.9) +
    geom_line(aes(group = Dataset), linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_size_continuous(name = "% Expressing", limits = c(0, 100), range = c(1, 6)) +
    labs(
      title = paste("Avg Astrocyte Expression Rank -", g),
      x = "Developmental Age",
      y = "Scaled Avg Rank (0–1)",
      color = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )

  # Save to output_dir
  ggsave(
    filename = file.path(output_dir, paste0("rank_plot_AST_", g, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
}

#### 5. all cell types ####

# Get all cell types present in datasets
all_b_types <- unique(b@meta.data$Cell_type_raw)
all_v_types <- unique(v@meta.data$Lineage)

all_b <- process_rank_and_percent(
  rank_mat = rank_mat_b,
  expr_mat = expr_mat_b,
  seurat_obj = b,
  age_col = "Developmental_week",
  cell_col = "Cell_type_raw",
  cell_types = all_b_types,
  dataset_label = "Braun"
)

# Process Velmeshev non-neuronal
all_v <- process_rank_and_percent(
  rank_mat = rank_mat_v,
  expr_mat = expr_mat_v,
  seurat_obj = v,
  age_col = "Age_Range",
  cell_col = "Lineage",
  cell_types = all_v_types,
  dataset_label = "Velmeshev"
)

# Combine datasets
all_data <- rbind(all_b, all_v)

# Ensure Age is an ordered factor (reuse age_levels vector)
all_data$Age <- factor(all_data$Age, levels = age_levels, ordered = TRUE)

# --- Plotting loop for non-neuronal ---
for (g in genes_of_interest) {
  gene_data <- all_data[all_data$Gene == g, ]
  
  # Aggregate by Gene + Age + Dataset
  agg_rank <- aggregate(avg_rank ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  agg_pct <- aggregate(percent_expressing ~ Gene + Age + Dataset, data = gene_data, FUN = mean, na.rm = TRUE)
  
  agg_df <- merge(agg_rank, agg_pct, by = c("Gene", "Age", "Dataset"))
  
  # Plot
  p <- ggplot(agg_df, aes(x = Age, y = avg_rank, color = Dataset, size = percent_expressing)) +
    geom_point(alpha = 0.9) +
    geom_line(aes(group = Dataset), linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_size_continuous(name = "% Expressing", limits = c(0, 100), range = c(1, 6)) +
    labs(
      title = paste("Avg All Cell Expression Rank -", g),
      x = "Developmental Age",
      y = "Scaled Avg Rank (0–1)",
      color = "Dataset"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save plot
  ggsave(
    filename = file.path(output_dir, paste0("rank_plot_all_", g, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
}

#### check for batch effects ####

# Check that the rownames match
all(rownames(rank_mat_b) == rownames(rank_mat_v))  # should return TRUE

# Combine them column-wise
combined_rank_matrix <- cbind(rank_mat_b, rank_mat_v)

  # Extract metadata from each Seurat object
meta_b <- b@meta.data
meta_v <- v@meta.data

# Add unique cell names if needed (to avoid duplicate column names)
colnames(rank_mat_b) <- paste0("cellb_", colnames(rank_mat_b))
colnames(rank_mat_v) <- paste0("cellv_", colnames(rank_mat_v))
rownames(meta_b) <- colnames(rank_mat_b)
rownames(meta_v) <- colnames(rank_mat_v)

# Combine metadata
meta_combined <- bind_rows(meta_b, meta_v)

# Create combined metadata columns
meta_combined$Developmental_Age <- paste(meta_combined$Developmental_week, meta_combined$Age_Range, sep = "_")
meta_combined$Cell_Type <- paste(meta_combined$Subclass, meta_combined$Lineage, sep = "_")
meta_combined$Region <- paste(meta_combined$Subdivision, meta_combined$Region_Broad, sep = "_")

# Create dummy count matrix (zeroes), same dimensions as your rank matrix
tmp <- matrix(0, nrow = nrow(combined_rank_matrix), ncol = ncol(combined_rank_matrix))
rownames(tmp) <- rownames(combined_rank_matrix)
colnames(tmp) <- colnames(combined_rank_matrix)

# Create Seurat object
seurat_combined <- CreateSeuratObject(counts = tmp, meta.data = meta_combined)

# Add the ranked data to the data slot (NOT counts)
seurat_combined <- SetAssayData(seurat_combined, assay = "RNA", slot = "data", new.data = combined_rank_matrix)

### generate UMAP ###
seurat_combined <- ScaleData(seurat_combined, verbose = FALSE)
seurat_combined <- RunPCA(seurat_combined, features = rownames(seurat_combined))

seurat_combined <- RunUMAP(seurat_combined, dims = 1:20)
umap_plot <- DimPlot(seurat_combined, reduction = "umap", group.by = "Developmental_Age")
ggsave("umap_plot_by_dataset.pdf", umap_plot, width = 6, height = 5)                          

saveRDS(seurat_combined, file = "seurat_combined_with_umap.rds")

                          

