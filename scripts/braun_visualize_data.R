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
library(readxl)

# #### 0. preprocessing 
# seurat_obj_path <- file.path(base_dir, "seurat_obj_merged_layers_only.rds")
# seurat_obj <- readRDS(seurat_obj_path)
# 
# # log transform fr this time
# expr_data <- GetAssayData(seurat_obj, assay = "merged", slot = "data")
# log_expr <- log2(expr_data + 1)
# log_expr_sparse <- as(log_expr, "CsparseMatrix")
# seurat_obj <- SetAssayData(seurat_obj, assay = "merged", slot = "data", new.data = log_expr_sparse)
# 
# # subset to 266 common asd risk genes
# common_asd_risk_genes <- readRDS("C:/Abbie/research/seurat/common_asd_risk_genes.rds")
# seurat_obj <- seurat_obj[common_asd_risk_genes, ]
# 
# # make dim reduction obj
# tsne_coords <- as.matrix(seurat_obj@meta.data[, c("tSNE_x", "tSNE_y")])
# colnames(tsne_coords) <- c("tSNE_1", "tSNE_2")
# 
# tsne_reduction <- CreateDimReducObject(
#   embeddings = tsne_coords,
#   key = "TSNE_",
#   assay = DefaultAssay(seurat_obj)
# )
# 
# seurat_obj@reductions$tsne <- tsne_reduction
# 
# 
# # sort age ranges from earliest to latest
# age_levels <- c("5", "5.5", "6", "6.6", "6.7", "6.9", "7", "7.5", "8", "8.1", "8.5", "9.2", "9.5", "10", "11.5", "12", "13", "14")
# seurat_obj$Developmental_week <- factor(seurat_obj$Developmental_week, levels = age_levels)
# 
# 
#  # add subclass labels for cluster annotations
# 
# cluster_anno <- read.csv(file.path(base_dir, "table_S2.csv"), header = TRUE)
# preac_clusters <- cluster_anno$PoolCleanOrder[grepl("PREAC", cluster_anno$AutoAnnotation)]
# meta <- seurat_obj@meta.data
# 
# cluster_anno$PoolCleanOrder <- as.character(cluster_anno$PoolCleanOrder)
# seurat_obj$Cell_clusters <- as.character(seurat_obj$Cell_clusters)
# cluster_to_subclass <- setNames(cluster_anno$Subclass, cluster_anno$PoolCleanOrder)
# subclass_vector <- cluster_to_subclass[seurat_obj$Cell_clusters]
# seurat_obj@meta.data$Subclass <- subclass_vector
# 
# seurat_obj$Subclass[seurat_obj$Cell_clusters %in% preac_clusters] <- "PREAC"
# 
# 
# saveRDS(seurat_obj, file.path(base_dir, "seurat_obj_merged_layers_only_processed.rds"))


seurat_obj_path <- file.path(base_dir, "seurat_obj_merged_layers_only_processed.rds")
seurat_obj <- readRDS(seurat_obj_path)

genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")
group_vars <- c("Lineage", "Region", "Subclass")

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


#### 2. Dot Plot ####

### 2a. Avg expression & % expressing ###
for (group_var in group_vars) {
  filename <- paste0("dot_plot_", tolower(group_var), ".pdf")
  filepath <- file.path(output_dir, filename)
  
  if (!file.exists(filepath)) {
    message("Generating and saving: ", filepath)
    
    title_text <- paste(
      paste(genes_of_interest, collapse = ", "),
      "expression by",
      tolower(gsub("_", " ", group_var))
    )
    
    p <- DotPlot(
      object = seurat_obj,
      features = genes_of_interest,
      group.by = group_var,
      cols = c("grey", "blue")
    ) +
      ggtitle(title_text) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    ggsave(filepath, plot = p, width = 10, height = 6, units = "in")
  } else {
    message("Skipped (already exists): ", filepath)
  }
}

### 2b. Avg rank & % expressing ###

for (group_var in group_vars) {
  
  # File paths
  file_rank <- file.path(output_dir, paste0("avg_rank_", tolower(group_var), ".xlsx"))
  file_pct_expr <- file.path(output_dir, paste0("pct_expr_", tolower(group_var), ".xlsx"))
  filename <- paste0("dot_plot_avg_rank_", tolower(group_var), ".pdf")
  filepath <- file.path(output_dir, filename)
  
  # Skip if file already exists
  if (!file.exists(filepath)) {
    message("Generating and saving: ", filepath)
    
    rank_df <- read_xlsx(file_rank) %>% as.data.frame()
    pct_df <- read_xlsx(file_pct_expr) %>% as.data.frame()
    rank_df$gene <- read_xlsx(file_rank)$gene
    pct_df$gene <- read_xlsx(file_pct_expr)$gene
    
    # Pivot to long format
    rank_long <- pivot_longer(rank_df, -gene, names_to = "celltype", values_to = "avg_rank")
    pct_long  <- pivot_longer(pct_df, -gene, names_to = "celltype", values_to = "pct_expr")
    
    dot_df <- left_join(rank_long, pct_long, by = c("gene", "celltype"))

    dot_df_subset <- dot_df %>% filter(gene %in% genes_of_interest)
    
    # # Handle missing/NA value
    # dot_df_subset$pct_expr[is.na(dot_df_subset$pct_expr)] <- 0
    # dot_df_subset$avg_rank[is.na(dot_df_subset$avg_rank)] <- max(dot_df_subset$avg_rank, na.rm = TRUE)
  
    # plot
    title_text <- paste(
      paste(genes_of_interest, collapse = ", "),
      "ranked expression by",
      gsub("_", " ", tolower(group_var))
    )
    
    p <- ggplot(dot_df_subset, aes(x = gene, y = celltype)) +
      geom_point(aes(size = pct_expr, color = avg_rank)) +
      scale_color_viridis_c(direction = -1, name = "Avg Rank") +
      scale_size(range = c(1, 8), name = "% Expressing") +
      ggtitle(title_text) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank()
      )
    
    # save
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    ggsave(filepath, plot = p, width = 10, height = 6, units = "in")
  } else {
    message("Skipped (already exists): ", filepath)
  }
}


#### 3. Scatter Plot ####

lineages <- c("All", "Erythrocyte", "Fibroblast", "Glioblast", "Immune", 
              "Neuroblast", "Neuron", "Neuronal IPC", "Oligo", "Radial glia", 
              "Vascular", "Neural crest", "Placodes")
age_levels <- c("5", "5.5", "6", "6.6", "6.7", "6.9", "7", "7.5", "8", "8.1", 
                "8.5", "9.2", "9.5", "10", "11.5", "12", "13", "14")
seurat_obj$Developmental_week <- factor(seurat_obj$Developmental_week, levels = age_levels)

# helper function to get summary data
get_summary_df <- function(seurat_obj, gene, lineage, age_levels) {
  if (lineage == "All") {
    seurat_sub <- seurat_obj
  } else {
    cells <- WhichCells(seurat_obj, expression = Cell_type_raw == lineage)
    seurat_sub <- subset(seurat_obj, cells = cells)
  }
  
  exprs_data <- FetchData(seurat_sub, vars = c(gene, "Developmental_week"))
  
  summary_df <- exprs_data %>%
    mutate(Developmental_week = factor(Developmental_week, levels = age_levels)) %>%
    group_by(Developmental_week) %>%
    summarize(
      avg_exprs = mean(.data[[gene]], na.rm = TRUE),
      percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100,
      .groups = "drop"
    ) %>%
    complete(Developmental_week = age_levels, fill = list(avg_exprs = 0, percent_exprs = 0)) %>%
    mutate(Gene = gene, Cell_type_raw = lineage)
  
  return(summary_df)
}

# prepare all summaries
all_summaries <- list()
for (gene in genes_of_interest) {
  for (lineage in lineages) {
    sum_df <- get_summary_df(seurat_obj, gene, lineage, age_levels)
    all_summaries[[paste(gene, lineage, sep = "_")]] <- sum_df
  }
}
combined_summary <- bind_rows(all_summaries)

# # find top 3 max avg. expressions
# top_expr_rows <- combined_summary %>%
#   arrange(desc(avg_exprs)) %>%
#   slice(1:3)
# 
# excluded_gene <- top_expr_rows$Gene[1]
# excluded_lineage <- top_expr_rows$Lineage[1]
# 
# # remove the top max row for consistent y-scale
# target_summary <- combined_summary %>%
#   filter(!(Gene == excluded_gene & Lineage == excluded_lineage))

# global axis scales
y_max <- ceiling(max(combined_summary$avg_exprs, na.rm = TRUE) * 10) / 10
max_percent <- ceiling(max(combined_summary$percent_exprs, na.rm = TRUE))

# function to generate & save plot if new
generate_and_save_age_expr_plot <- function(df, gene, lineage, y_max, max_percent, output_dir) {
  lineage_clean <- tolower(gsub("[^a-zA-Z0-9]", "", lineage))
  gene_clean <- gene  # preserve case
  
  filename <- paste0("age_expr_", gene_clean, "_", lineage_clean, ".pdf")
  filepath <- file.path(output_dir, filename)
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  if (!file.exists(filepath)) {
    message("Generating and saving: ", filepath)
    
    subtitle_text <- ifelse(lineage == "All", "All Cells", paste(lineage, "Cells"))
    
    df <- df %>% mutate(Developmental_week = factor(Developmental_week, levels = age_levels))
    
    p <- ggplot(df, aes(x = Developmental_week, y = avg_exprs, group = 1)) +
      geom_line(color = "#2c7fb8", linewidth = 1) +
      geom_point(aes(size = percent_exprs), color = "#2c7fb8") +
      scale_size_continuous(limits = c(0, max_percent), breaks = scales::pretty_breaks(n = 5)) +
      theme_minimal() +
      labs(
        title = paste("Expression of", gene, "across Developmental Weeks"),
        subtitle = subtitle_text,
        x = "Developmental Week",
        y = "Average Expression",
        size = "% Expressing Cells"
      ) +
      ylim(0, y_max) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(filepath, plot = p, width = 8, height = 5, units = "in")
  } else {
    message("Skipped (already exists): ", filepath)
  }
}

# loop to generate plots
for (gene in genes_of_interest) {
  for (lineage in lineages) {
  
    df <- combined_summary %>%
      filter(Gene == gene, Cell_type_raw == lineage)
    
    generate_and_save_age_expr_plot(df, gene, lineage, y_max, max_percent, output_dir)
  }
}

#### 4. dataset statistics ####

# General save function: generates and saves if file doesn't exist
generate_and_save_data_stats_plot <- function(plot_obj, filename, width = 8, height = 5) {
  filepath <- file.path(output_dir, filename)
  if (!file.exists(filepath)) {
    message("Saving: ", filepath)
    ggsave(filepath, plot_obj, width = width, height = height)
  } else {
    message("Skipped (already exists): ", filepath)
  }
}

# 4a. cell counts per donor
donor_counts <- as.data.frame(table(seurat_obj$Individual))
colnames(donor_counts) <- c("Individual", "Cell_Count")

p1 <- ggplot(donor_counts, aes(x = Individual, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Individual", y = "Number of Cells", title = "Cell Counts per Donor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

generate_and_save_data_stats_plot(p1, "cell_counts_per_donor.pdf", width = 10)

### 4b. cell counts per age range ###
age_counts <- as.data.frame(table(seurat_obj$Developmental_week))
colnames(age_counts) <- c("Developmental_week", "Cell_Count")

p2 <- ggplot(age_counts, aes(x = Developmental_week, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Developmental Week", y = "Number of Cells", title = "Developmental Week") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p2, "cell_counts_per_developmental_week.pdf")

### 4c. cell counts per region ###
region_counts <- as.data.frame(table(seurat_obj$Region))
colnames(region_counts) <- c("Region", "Cell_Count")

p3 <- ggplot(region_counts, aes(x = Region, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Region", y = "Number of Cells", title = "Cell Counts per Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p3, "cell_counts_per_region.pdf")

### 4d. cell counts per region -- specific ###
subdiv_counts <- as.data.frame(table(seurat_obj$Subdivision))
colnames(subdiv_counts) <- c("Subdivision", "Cell_Count")

p4 <- ggplot(subdiv_counts, aes(x = Subdivision, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Subdivision", y = "Number of Cells", title = "Cell Counts per Subdivision") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p4, "cell_counts_per_subdivision.pdf")

### 4e. cell counts per lineage ###
lineage_counts <- as.data.frame(table(seurat_obj$Cell_type_raw))
colnames(lineage_counts) <- c("Cell_type_raw", "Cell_Count")

p5 <- ggplot(lineage_counts, aes(x = Cell_type_raw, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Cell_type_raw", y = "Number of Cells", title = "Cell Counts per Cell Type Raw") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p5, "cell_counts_per_cell_type_raw.pdf")

### 4f. cell type composition per region (%) ###
region_celltype_counts <- seurat_obj@meta.data %>%
  group_by(Region, Cell_type_raw) %>%
  summarise(Count = n(), .groups = "drop")

region_celltype_percent <- region_celltype_counts %>%
  group_by(Region) %>%
  mutate(Percent = Count / sum(Count) * 100)

p6 <- ggplot(region_celltype_percent, aes(x = Region, y = Percent, fill = Cell_type_raw)) +
  geom_bar(stat = "identity") +
  labs(x = "Region", y = "Cell Type Composition (%)", fill = "Cell Type",
       title = "Cell Type Composition per Brain Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p6, "cell_type_composition_per_region.pdf", width = 10)



#### 5. Heatmap Analysis of New and Known Autism Risk Genes ####


# Load expression data and metadata
exprs_data_asd <- GetAssayData(seurat_obj, assay = "merged", slot = "data")
meta <- seurat_obj@meta.data
meta$cell <- rownames(meta)

# Load rank data
rank_data <- readRDS(file.path(base_dir, "cell_rankings_braun.rds"))

# Load ASD risk genes
asd_risk_genes <- read.csv("autism_risk_genes_combined.csv", sep = ",", header = TRUE)
asd_risk_genes <- asd_risk_genes[asd_risk_genes$Gene %in% rownames(seurat_obj), ]
gene_status <- setNames(asd_risk_genes$Status, asd_risk_genes$Gene)

# Set output directory and file paths
file_exprs_region <- file.path(output_dir, "avg_expr_region.xlsx")
file_exprs_lineage <- file.path(output_dir, "avg_expr_lineage.xlsx")
file_exprs_subclass <- file.path(output_dir, "avg_expr_subclass.xlsx")
file_heatmap_region <- file.path(output_dir, "heatmap_avg_expr_region.pdf")
file_heatmap_lineage <- file.path(output_dir, "heatmap_avg_expr_lineage.pdf")
file_heatmap_subclass <- file.path(output_dir, "heatmap_avg_expr_subclass.pdf")
file_pct_exprs_region <- file.path(output_dir, "pct_expr_region.xlsx")
file_pct_exprs_lineage <- file.path(output_dir, "pct_expr_lineage.xlsx")
file_pct_exprs_subclass <- file.path(output_dir, "pct_expr_subclass.xlsx")
file_heatmap_pct_region <- file.path(output_dir, "heatmap_pct_expr_region.pdf")
file_heatmap_pct_lineage <- file.path(output_dir, "heatmap_pct_expr_lineage.pdf")
file_heatmap_pct_subclass <- file.path(output_dir, "heatmap_pct_expr_subclass.pdf")
file_rank_region <- file.path(output_dir, "avg_rank_region.xlsx")
file_rank_lineage <- file.path(output_dir, "avg_rank_lineage.xlsx")
file_rank_subclass <- file.path(output_dir, "avg_rank_subclass.xlsx")
file_heatmap_rank_region <- file.path(output_dir, "heatmap_avg_rank_region.pdf")
file_heatmap_rank_lineage <- file.path(output_dir, "heatmap_avg_rank_lineage.pdf")
file_heatmap_rank_subclass <- file.path(output_dir, "heatmap_avg_rank_subclass.pdf")


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
gc()

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
gc()

### 5c. avg expression per gene by subclass ###
if (!file.exists(file_exprs_subclass)) {
  message("Generating average expression table by Subclass")
  
  df_subclass <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Subclass")], by = "cell")
  
  avg_exprs_subclass <- df_subclass %>%
    group_by(gene, Subclass) %>%
    summarise(avg_exprs = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = Subclass, values_from = avg_exprs)
  
  write_xlsx(avg_exprs_subclass, file_exprs_subclass)
} else {
  message("Average expression table by Subclass exists, loading...")
  avg_exprs_subclass <- readxl::read_xlsx(file_exprs_subclass)
}
gc()

### 5d. % expression per gene by region ###
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
gc()

### 5e. % expression per gene by lineage ###
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
gc()

### 5f. % expression per gene by subclass ###
if (!file.exists(file_pct_exprs_subclass)) {
  message("Generating % expression table by Subclass")
  
  df_subclass <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Subclass")], by = "cell")
  
  pct_exprs_subclass <- df_subclass %>%
    group_by(gene, Subclass) %>%
    summarise(pct_exprs = sum(expression > 0) / n() * 100, .groups = "drop") %>%
    pivot_wider(names_from = Subclass, values_from = pct_exprs)
  
  write_xlsx(pct_exprs_subclass, file_pct_exprs_subclass)
} else {
  message("% Expression table by Subclass exists, loading...")
  pct_exprs_subclass <- readxl::read_xlsx(file_pct_exprs_subclass)
}
gc()

### 5g. avg rank per gene by region ###

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
gc()


### 5h. avg rank per gene by lineage ###

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
gc()

### 5i. avg rank per gene by subclass ###

if (!file.exists(file_rank_subclass)) {
  message("Generating average rank table by Subclass")
  
  df_subclass <- rank_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "rank") %>%
    left_join(meta[, c("cell", "Subclass")], by = "cell")
  
  avg_rank_subclass <- df_subclass %>%
    group_by(gene, Subclass) %>%
    summarise(avg_rank = mean(rank), .groups = "drop") %>%
    pivot_wider(names_from = Subclass, values_from = avg_rank)
  
  write_xlsx(avg_rank_subclass, file_rank_subclass)
} else {
  message("Average rank table by Subclass exists, loading...")
  avg_rank_subclass <- readxl::read_xlsx(file_rank_subclass)
}
gc()

### 5j. heatmap by lineage ###
if (!file.exists(file_heatmap_lineage)) {
  message("Generating heatmap by Cell_type_raw")
  
  exprs_mat_lineage <- as.matrix(avg_exprs_lineage[, -1])
  rownames(exprs_mat_lineage) <- avg_exprs_lineage$gene
  exprs_mat_lineage_t <- t(exprs_mat_lineage)
  
  gene_status <- gene_status[colnames(exprs_mat_lineage_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
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
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_lineage, width = 20, height = 8)
  draw(ht_lineage)
  dev.off()
} else {
  message("Cell_type_raw heatmap file exists, skipping generation.")
}

### 5k. heatmap by region ###
if (!file.exists(file_heatmap_region)) {
  message("Generating heatmap by Region")
  
  exprs_mat_region <- as.matrix(avg_exprs_region[, -1])
  rownames(exprs_mat_region) <- avg_exprs_region$gene
  exprs_mat_region_t <- t(exprs_mat_region)
  
  gene_status <- gene_status[colnames(exprs_mat_region_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
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
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_region, width = 20, height = 8)
  draw(ht_region)
  dev.off()
} else {
  message("Region heatmap file exists, skipping generation.")
}

### 5l. heatmap by subclass ###
if (!file.exists(file_heatmap_subclass)) {
  message("Generating heatmap by Subclass")
  
  exprs_mat_subclass <- as.matrix(avg_exprs_subclass[, -1])
  rownames(exprs_mat_subclass) <- avg_exprs_subclass$gene
  exprs_mat_subclass_t <- t(exprs_mat_subclass)
  
  gene_status <- gene_status[colnames(exprs_mat_subclass_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
  exprs_colors <- colorRamp2(
    c(min(exprs_mat_subclass_t), median(exprs_mat_subclass_t), max(exprs_mat_subclass_t)),
    c("blue", "white", "red")
  )
  
  ht_subclass <- Heatmap(
    exprs_mat_subclass_t,
    name = "Avg Expression",
    column_title = "Gene",
    row_title = "Subclass",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_subclass, width = 20, height = 8)
  draw(ht_subclass)
  dev.off()
} else {
  message("Subclass heatmap file exists, skipping generation.")
}

### 5m. heatmap of % exprs by lineage ###
if (!file.exists(file_heatmap_pct_lineage)) {
  message("Generating % exprs heatmap by Cell_type_raw")
  
  pct_mat_lineage <- as.matrix(pct_exprs_lineage[, -1])
  rownames(pct_mat_lineage) <- pct_exprs_lineage$gene
  pct_mat_lineage_t <- t(pct_mat_lineage)
  
  gene_status <- gene_status[colnames(pct_mat_lineage_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
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
    col = pct_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_pct_lineage, width = 20, height = 8)
  draw(ht_pct_lineage)
  dev.off()
} else {
  message("% Exprs Cell_type_raw heatmap exists, skipping generation.")
}

### 5n. heatmap of % exprs by region ###
if (!file.exists(file_heatmap_pct_region)) {
  message("Generating % exprs heatmap by Region")
  
  pct_mat_region <- as.matrix(pct_exprs_region[, -1])
  rownames(pct_mat_region) <- pct_exprs_region$gene
  pct_mat_region_t <- t(pct_mat_region)
  
  gene_status <- gene_status[colnames(pct_mat_region_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
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
    col = pct_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_pct_region, width = 20, height = 8)
  draw(ht_pct_region)
  dev.off()
} else {
  message("% Exprs Region heatmap exists, skipping generation.")
}

### 5o. heatmap of % exprs by subclass ###
if (!file.exists(file_heatmap_pct_subclass)) {
  message("Generating % exprs heatmap by Subclass")
  
  pct_mat_subclass <- as.matrix(pct_exprs_subclass[, -1])
  rownames(pct_mat_subclass) <- pct_exprs_subclass$gene
  pct_mat_subclass_t <- t(pct_mat_subclass)
  
  gene_status <- gene_status[colnames(pct_mat_subclass_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
  pct_colors <- colorRamp2(
    c(min(pct_mat_subclass_t), median(pct_mat_subclass_t), max(pct_mat_subclass_t)),
    c("blue", "white", "red")
  )
  
  ht_pct_subclass <- Heatmap(
    pct_mat_subclass_t,
    name = "% Expressing",
    column_title = "Gene",
    row_title = "Subclass",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = pct_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_pct_subclass, width = 20, height = 8)
  draw(ht_pct_subclass)
  dev.off()
} else {
  message("% Exprs Subclass heatmap exists, skipping generation.")
}

### 5p. heatmap for rank by region ###
if (!file.exists(file_heatmap_rank_region)) {
  message("Generating rank heatmap by Region")
  
  rank_mat_region <- as.matrix(avg_rank_region[, -1])
  rownames(rank_mat_region) <- avg_rank_region$gene
  rank_mat_region_t <- t(rank_mat_region)
  
  gene_status <- gene_status[colnames(rank_mat_region_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
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
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_rank_region, width = 20, height = 8)
  draw(ht_region)
  dev.off()
} else {
  message("Region rank heatmap file exists, skipping generation.")
}

### 5q. heatmap for rank by lineage ###
if (!file.exists(file_heatmap_rank_lineage)) {
  message("Generating rank heatmap by Cell_type_raw")
  
  rank_mat_lineage <- as.matrix(avg_rank_lineage[, -1])
  rownames(rank_mat_lineage) <- avg_rank_lineage$gene
  rank_mat_lineage_t <- t(rank_mat_lineage)
  
  gene_status <- gene_status[colnames(rank_mat_lineage_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
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
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_rank_lineage, width = 20, height = 8)
  draw(ht_lineage)
  dev.off()
} else {
  message("Cell_type_raw rank heatmap file exists, skipping generation.")
}

### 5r. heatmap for rank by subclass ###
if (!file.exists(file_heatmap_rank_subclass)) {
  message("Generating rank heatmap by Subclass")
  
  rank_mat_subclass <- as.matrix(avg_rank_subclass[, -1])
  rownames(rank_mat_subclass) <- avg_rank_subclass$gene
  rank_mat_subclass_t <- t(rank_mat_subclass)
  
  gene_status <- gene_status[colnames(rank_mat_subclass_t)] # reorder gene status
  # col annotation
  col_annot <- HeatmapAnnotation(
    GeneStatus = gene_status,
    col = list(GeneStatus = c("Known" = "steelblue", "New" = "salmon")),
    annotation_name_side = "left"
  )
  
  exprs_colors <- colorRamp2(
    c(min(rank_mat_subclass_t), median(rank_mat_subclass_t), max(rank_mat_subclass_t)),
    c("blue", "white", "red")
  )
  
  ht_subclass <- Heatmap(
    rank_mat_subclass_t,
    name = "Avg Rank",
    column_title = "Gene",
    row_title = "Subclass",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 6),
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_rank_subclass, width = 20, height = 8)
  draw(ht_subclass)
  dev.off()
} else {
  message("Subclass rank heatmap file exists, skipping generation.")
}

