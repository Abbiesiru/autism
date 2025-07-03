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

genes_of_interest <- c("SORCS1", "SORCS3")
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


