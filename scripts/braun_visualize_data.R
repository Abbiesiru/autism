base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/braun",
                   "gauss" = "/home/abbiew/single_cell/braun",
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
library(Seurat)
library(patchwork)
library(Matrix)
library(ggplot2)
library(data.table)
library(writexl)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(grid)

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_combined.rds")
seurat_obj <- readRDS(seurat_obj_path)

genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")
group_vars <- c("Cell_type_raw", "Region")

#### 0. make dim reduction object & sort dev weeks ####

tsne_coords <- as.matrix(seurat_obj@meta.data[, c("tSNE_x", "tSNE_y")])
colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")  

tsne_reduction <- CreateDimReducObject(
  embeddings = tsne_coords,
  key = "TSNE_",
  assay = DefaultAssay(seurat_obj)
)

seurat_obj@reductions$tsne <- tsne_reduction

# DimPlot(seurat_obj, reduction = "tsne", group.by = "Region") 


# sort age ranges from earliest to latest
age_levels <- c("5", "5.5", "6", "6.6", "6.7", "6.9", "7", "7.5", "8", "8.1", "8.5", "9.2", "9.5", "10", "11.5", "12", "13", "14")
seurat_obj$Developmental_week <- factor(seurat_obj$Developmental_week, levels = age_levels)

#### 1. Feature Plot ####

### 1a. create feature plot with annotations ###
generate_tsne_plot <- function(seurat_obj, group_var, gene, split_by_age = TRUE) {
  # Calculate centroids (for non-split case)
  tsne_df <- data.frame(
    x = seurat_obj@reductions$tsne@cell.embeddings[, 1],
    y = seurat_obj@reductions$tsne@cell.embeddings[, 2],
    Group = seurat_obj[[group_var]][, 1]
  )
  
  centroids <- tsne_df %>%
    group_by(Group) %>%
    summarize(x = mean(x), y = mean(y), .groups = "drop")
  
  # Determine consistent scale across plots
  expr_vals <- FetchData(seurat_obj, vars = gene)[, 1]
  min_expr <- min(expr_vals)
  max_expr <- max(expr_vals)
  
  if (split_by_age) {
    age_levels <- unique(seurat_obj$Developmental_week)
    
    plots <- lapply(age_levels, function(age) {
      subset_obj <- subset(seurat_obj, subset = Developmental_week == age)
      FeaturePlot(
        subset_obj,
        reduction = "tsne",
        features = gene,
        min.cutoff = min_expr,
        max.cutoff = max_expr
      ) +
        ggtitle(paste0("Week ", age))
      +
        scale_color_gradientn(
          colors = c("lightgrey", "blue"),
          limits = c(min_expr, max_expr),
          oob = scales::squish
        )
    })
    
    return(wrap_plots(plots, ncol = 2) +
             patchwork::plot_annotation(title = paste(group_var, "-", gene)))
    
  } else {
    FeaturePlot(
      seurat_obj,
      reduction = "tsne",
      features = gene,
      min.cutoff = min_expr,
      max.cutoff = max_expr
    ) +
      geom_text(
        data = centroids,
        aes(x = x, y = y, label = Group),
        size = 3,
        color = "black",
        inherit.aes = FALSE
      ) +
      ggtitle(paste("All Ages -", group_var, "-", gene))
    +
      scale_color_gradientn(
        colors = c("lightgrey", "blue"),
        limits = c(min_expr, max_expr),
        oob = scales::squish
      )
  }
}



### 1b. generate and save plot only if not already saved ###
generate_and_save_if_new <- function(seurat_obj, group_var, gene, split_by_age, folder, filename, base_width = 16, base_height = 6) {
  filepath <- file.path(folder, filename)
  
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  if (!file.exists(filepath)) {
    message("Generating and saving: ", filepath)
    
    if (split_by_age) {
      # count number of age levels
      n_panels <- length(unique(seurat_obj$Developmental_week))
      ncol <- 2
      nrow <- ceiling(n_panels / ncol)
      height <- base_height * nrow  # height grows with number of rows
    } else {
      height <- base_height * 2  # fixed for 1 plot
    }
    
    p <- generate_tsne_plot(seurat_obj, group_var = group_var, gene = gene, split_by_age = split_by_age)
    ggsave(filepath, plot = p, width = base_width, height = height, units = "in")
  } else {
    message("Skipped (already exists): ", filepath)
  }
}



### generate 8 plots: 2 genes x 2 annotations, split/unsplit ###

group_vars <- c("Cell_type_raw", "Region")
split_options <- c(TRUE, FALSE)

for (group_var in group_vars) {
  for (gene in genes_of_interest) {
    for (split in split_options) {
      split_label <- if (split) "split" else "all"
      filename <- paste0("feature_plot_", gene, "_", tolower(group_var), "_", split_label, ".pdf")
      
      generate_and_save_if_new(
        seurat_obj = seurat_obj,
        group_var = group_var,
        gene = gene,
        split_by_age = split,
        folder = output_dir,
        filename = filename,
      )
    }
  }
}


#### 5. Heatmap Analysis of New and Known Autism Risk Genes ####

# Load expression data and metadata
exprs_data_asd <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
meta <- seurat_obj@meta.data
meta$cell <- rownames(meta)

# Load rank data
rank_data <- readRDS(file.path(base_dir, "cell_rankings_braun.rds"))

# Set output directory and file paths
file_exprs_region <- file.path(output_dir, "avg_expr_region.xlsx")
file_exprs_lineage <- file.path(output_dir, "avg_expr_lineage.xlsx")
file_heatmap_region <- file.path(output_dir, "heatmap_avg_expr_region.pdf")
file_heatmap_lineage <- file.path(output_dir, "heatmap_avg_expr_lineage.pdf")
file_pct_exprs_region <- file.path(output_dir, "pct_expr_region.xlsx")
file_pct_exprs_lineage <- file.path(output_dir, "pct_expr_lineage.xlsx")
file_heatmap_pct_region <- file.path(output_dir, "heatmap_pct_expr_region.pdf")
file_heatmap_pct_lineage <- file.path(output_dir, "heatmap_pct_expr_lineage.pdf")
file_rank_region <- file.path(output_dir, "avg_rank_region.xlsx")
file_rank_lineage <- file.path(output_dir, "avg_rank_lineage.xlsx")
file_heatmap_rank_region <- file.path(output_dir, "heatmap_avg_rank_region.pdf")
file_heatmap_rank_lineage <- file.path(output_dir, "heatmap_avg_rank_lineage.pdf")

### 5a. avg expression per gene by region ###
if (!file.exists(file_exprs_region)) {
  message("Generating average expression table by Region")
  
  df_region <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Region")], by = "cell")
  
  avg_exprs_region <- df_region %>%
    group_by(gene, Region) %>%
    summarise(avg_exprs = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = Region, values_from = avg_exprs)
  
  write_xlsx(avg_exprs_region, file_exprs_region)
} else {
  message("Average expression table by Region exists, loading...")
  avg_exprs_region <- readxl::read_xlsx(file_exprs_region)
}

### 5b. avg expression per gene by Cell_type_raw ###
if (!file.exists(file_exprs_lineage)) {
  message("Generating average expression table by Cell_type_raw")
  
  df_lineage <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Cell_type_raw")], by = "cell")
  
  avg_exprs_lineage <- df_lineage %>%
    group_by(gene, Cell_type_raw) %>%
    summarise(avg_exprs = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = Cell_type_raw, values_from = avg_exprs)
  
  write_xlsx(avg_exprs_lineage, file_exprs_lineage)
} else {
  message("Average expression table by Cell_type_raw exists, loading...")
  avg_exprs_lineage <- readxl::read_xlsx(file_exprs_lineage)
}

### 5c. % expression per gene by region ###
if (!file.exists(file_pct_exprs_region)) {
  message("Generating % expression table by Region")
  
  df_region <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Region")], by = "cell")
  
  pct_exprs_region <- df_region %>%
    group_by(gene, Region) %>%
    summarise(pct_exprs = sum(expression > 0) / n() * 100, .groups = "drop") %>%
    pivot_wider(names_from = Region, values_from = pct_exprs)
  
  write_xlsx(pct_exprs_region, file_pct_exprs_region)
} else {
  message("% Expression table by Region exists, loading...")
  pct_exprs_region <- readxl::read_xlsx(file_pct_exprs_region)
}

### 5d. % expression per gene by lineage ###
if (!file.exists(file_pct_exprs_lineage)) {
  message("Generating % expression table by Cell_type_raw")
  
  df_lineage <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Cell_type_raw")], by = "cell")
  
  pct_exprs_lineage <- df_lineage %>%
    group_by(gene, Cell_type_raw) %>%
    summarise(pct_exprs = sum(expression > 0) / n() * 100, .groups = "drop") %>%
    pivot_wider(names_from = Cell_type_raw, values_from = pct_exprs)
  
  write_xlsx(pct_exprs_lineage, file_pct_exprs_lineage)
} else {
  message("% Expression table by Cell_type_raw exists, loading...")
  pct_exprs_lineage <- readxl::read_xlsx(file_pct_exprs_lineage)
}

### 5e. avg rank per gene by region ###

if (!file.exists(file_rank_region)) {
  message("Generating average rank table by Region")
  
  df_region <- rank_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "rank") %>%
    left_join(meta[, c("cell", "Region")], by = "cell")
  
  avg_rank_region <- df_region %>%
    group_by(gene, Region) %>%
    summarise(avg_rank = mean(rank), .groups = "drop") %>%
    pivot_wider(names_from = Region, values_from = avg_rank)
  
  write_xlsx(avg_rank_region, file_rank_region)
} else {
  message("Average rank table by Region exists, loading...")
  avg_rank_region <- readxl::read_xlsx(file_rank_region)
}

### 5f. avg rank per gene by lineage ###

if (!file.exists(file_rank_lineage)) {
  message("Generating average rank table by Cell_type_raw")
  
  df_lineage <- rank_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "rank") %>%
    left_join(meta[, c("cell", "Cell_type_raw")], by = "cell")
  
  avg_rank_lineage <- df_lineage %>%
    group_by(gene, Cell_type_raw) %>%
    summarise(avg_rank = mean(rank), .groups = "drop") %>%
    pivot_wider(names_from = Cell_type_raw, values_from = avg_rank)
  
  write_xlsx(avg_rank_lineage, file_rank_lineage)
} else {
  message("Average rank table by Cell_type_raw exists, loading...")
  avg_rank_lineage <- readxl::read_xlsx(file_rank_lineage)
}

### 5g. heatmap by lineage ###
if (!file.exists(file_heatmap_lineage)) {
  message("Generating heatmap by Cell_type_raw")
  
  exprs_mat_lineage <- as.matrix(avg_exprs_lineage[, -1])
  rownames(exprs_mat_lineage) <- avg_exprs_lineage$gene
  exprs_mat_lineage_t <- t(exprs_mat_lineage)
  
  exprs_colors <- colorRamp2(
    c(min(exprs_mat_lineage_t), median(exprs_mat_lineage_t), max(exprs_mat_lineage_t)),
    c("blue", "white", "red")
  )
  
  ht_lineage <- Heatmap(
    exprs_mat_lineage_t,
    name = "Avg Expression",
    column_title = "Gene",
    row_title = "Cell Type",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = exprs_colors
  )
  
  pdf(file_heatmap_lineage, width = 20, height = 8)
  draw(ht_lineage)
  dev.off()
} else {
  message("Cell_type_raw heatmap file exists, skipping generation.")
}

### 5h. heatmap by region ###
if (!file.exists(file_heatmap_region)) {
  message("Generating heatmap by Region")
  
  exprs_mat_region <- as.matrix(avg_exprs_region[, -1])
  rownames(exprs_mat_region) <- avg_exprs_region$gene
  exprs_mat_region_t <- t(exprs_mat_region)
  
  exprs_colors <- colorRamp2(
    c(min(exprs_mat_region_t), median(exprs_mat_region_t), max(exprs_mat_region_t)),
    c("blue", "white", "red")
  )
  
  ht_region <- Heatmap(
    exprs_mat_region_t,
    name = "Avg Expression",
    column_title = "Gene",
    row_title = "Region",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = exprs_colors
  )
  
  pdf(file_heatmap_region, width = 20, height = 8)
  draw(ht_region)
  dev.off()
} else {
  message("Region heatmap file exists, skipping generation.")
}

### 5i. heatmap of % exprs by lineage ###
if (!file.exists(file_heatmap_pct_lineage)) {
  message("Generating % exprs heatmap by Cell_type_raw")
  
  pct_mat_lineage <- as.matrix(pct_exprs_lineage[, -1])
  rownames(pct_mat_lineage) <- pct_exprs_lineage$gene
  pct_mat_lineage_t <- t(pct_mat_lineage)
  
  pct_colors <- colorRamp2(
    c(min(pct_mat_lineage_t), median(pct_mat_lineage_t), max(pct_mat_lineage_t)),
    c("blue", "white", "red")
  )
  
  ht_pct_lineage <- Heatmap(
    pct_mat_lineage_t,
    name = "% Expressing",
    column_title = "Gene",
    row_title = "Cell Type",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = pct_colors
  )
  
  pdf(file_heatmap_pct_lineage, width = 20, height = 8)
  draw(ht_pct_lineage)
  dev.off()
} else {
  message("% Exprs Cell_type_raw heatmap exists, skipping generation.")
}

### 5j. heatmap of % exprs by region ###
if (!file.exists(file_heatmap_pct_region)) {
  message("Generating % exprs heatmap by Region")
  
  pct_mat_region <- as.matrix(pct_exprs_region[, -1])
  rownames(pct_mat_region) <- pct_exprs_region$gene
  pct_mat_region_t <- t(pct_mat_region)
  
  pct_colors <- colorRamp2(
    c(min(pct_mat_region_t), median(pct_mat_region_t), max(pct_mat_region_t)),
    c("blue", "white", "red")
  )
  
  ht_pct_region <- Heatmap(
    pct_mat_region_t,
    name = "% Expressing",
    column_title = "Gene",
    row_title = "Region",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = pct_colors
  )
  
  pdf(file_heatmap_pct_region, width = 20, height = 8)
  draw(ht_pct_region)
  dev.off()
} else {
  message("% Exprs Region heatmap exists, skipping generation.")
}


### 5k. heatmap for rank by region ###
if (!file.exists(file_heatmap_rank_region)) {
  message("Generating rank heatmap by Region")
  
  rank_mat_region <- as.matrix(avg_rank_region[, -1])
  rownames(rank_mat_region) <- avg_rank_region$gene
  rank_mat_region_t <- t(rank_mat_region)
  
  exprs_colors <- colorRamp2(
    c(min(rank_mat_region_t), median(rank_mat_region_t), max(rank_mat_region_t)),
    c("blue", "white", "red")
  )
  
  ht_region <- Heatmap(
    rank_mat_region_t,
    name = "Avg Rank",
    column_title = "Gene",
    row_title = "Region",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = exprs_colors
  )
  
  pdf(file_heatmap_rank_region, width = 20, height = 8)
  draw(ht_region)
  dev.off()
} else {
  message("Region rank heatmap file exists, skipping generation.")
}

### 5l. heatmap for rank by lineage ###
if (!file.exists(file_heatmap_rank_lineage)) {
  message("Generating rank heatmap by Cell_type_raw")
  
  rank_mat_lineage <- as.matrix(avg_rank_lineage[, -1])
  rownames(rank_mat_lineage) <- avg_rank_lineage$gene
  rank_mat_lineage_t <- t(rank_mat_lineage)
  
  exprs_colors <- colorRamp2(
    c(min(rank_mat_lineage_t), median(rank_mat_lineage_t), max(rank_mat_lineage_t)),
    c("blue", "white", "red")
  )
  
  ht_lineage <- Heatmap(
    rank_mat_lineage_t,
    name = "Avg Rank",
    column_title = "Gene",
    row_title = "Cell Type",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = exprs_colors
  )
  
  pdf(file_heatmap_rank_lineage, width = 20, height = 8)
  draw(ht_lineage)
  dev.off()
} else {
  message("Cell_type_raw rank heatmap file exists, skipping generation.")
}


