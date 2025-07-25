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


seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_common_genes.rds")
seurat_obj <- readRDS(seurat_obj_path)
genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")
group_vars <- c("Lineage")

#### 0. UMAP ####
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

brain_cell_types <- c("Astrocytes", 
                      "Excitatory neurons", 
                      "Glial progenitors", 
                      "Inhibitory neurons", 
                      "Microglia", 
                      "Oligodendrocytes", 
                      "Oligodendrocyte precursors" 
                      )
seurat_obj_brain <- subset(seurat_obj, subset = Lineage %in% brain_cell_types)

seurat_obj$Region_Broad <- recode(
  seurat_obj$Region_Broad,
  "CC" = "Cingulate cortex",
  "CGE" = "Caudal ganglionic eminence",
  "CTX" = "Cortex",
  "FC" = "Frontal/prefrontal cortex",
  "GE" = "Ganglionic eminence",
  "IC" = "Insular cortex",
  "LGE" = "Lateral ganglionic eminence",
  "MC" = "Motor cortex",
  "MGE" = "medial ganglionic eminence",
  "TC" = "Temporal cortex"
)


p <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "Age_Range",     # color by cell type
  label = TRUE,               # add labels to clusters
  label.size = 3,
  pt.size = 1.5,              # control point size
) + ggtitle("UMAP - Developmental Age")

ggsave(file.path(output_dir, "umap_plot_age.png"), p, width = 8, height = 6)

#### 1. Feature Plot ####

### 1a. create feature plot with annotations ###
age_order <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")

# Set Age_Range as factor with correct order
seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_order)

generate_umap_plot <- function(seurat_obj, group_var, gene, split_by_age = TRUE) {
  # Calculate centroids
  umap_df <- data.frame(
    x = seurat_obj@reductions$umap@cell.embeddings[, 1],
    y = seurat_obj@reductions$umap@cell.embeddings[, 2],
    Group = seurat_obj[[group_var]][, 1]
  )

  centroids <- umap_df %>%
    group_by(Group) %>%
    summarize(x = mean(x), y = mean(y), .groups = "drop")
  
  
  # If splitting by age
  if (split_by_age) {
    age_levels <- levels(seurat_obj$Age_Range)
    
    plots <- lapply(age_levels, function(age) {
      subset_obj <- subset(seurat_obj, subset = Age_Range == age)
      
      FeaturePlot(subset_obj, reduction = "umap", features = gene) +
        geom_text(
          data = centroids,
          aes(x = x, y = y, label = Group),
          size = 3,
          color = "black",
          inherit.aes = FALSE
        ) +
        ggtitle(as.character(age))
    })
    
    return(wrap_plots(plots, ncol = 4) + 
             patchwork::plot_annotation(title = paste(group_var, "-", gene)))  
    
  } else {
    # Single UMAP, not split by age
    FeaturePlot(seurat_obj, reduction = "umap", features = gene) +
      geom_text(
        data = centroids,
        aes(x = x, y = y, label = Group),
        size = 8,
        color = "black",
        inherit.aes = FALSE
      ) +
      ggtitle(paste("All Ages -", group_var, "-", gene))
  }
}

### 1b. generate and save plot only if not already saved ###
generate_and_save_if_new <- function(seurat_obj, group_var, gene, split_by_age, folder, filename, width = 16, height = 12) {
  # Construct full file path
  filepath <- file.path(folder, filename)
  
  # Create output folder if needed
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  # Check existence
  if (!file.exists(filepath)) {
    message("Generating and saving: ", filepath)
    p <- generate_umap_plot(seurat_obj, group_var = group_var, gene = gene, split_by_age = split_by_age)
    ggsave(filepath, plot = p, width = width, height = height, units = "in")
  } else {
    message("Skipped (already exists): ", filepath)
  }
}


### generate 12 plots: 3 genes x 2 annotations, split/unsplit ###

split_options <- c(FALSE)

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
        width = 16,
        height = 12
      )
    }
  }
}

#### 2. Dot Plot ####

sig_df <- fread(file.path(base_dir, "optimized_chisq_results_lineage.csv"))

### 2a. Avg expression & % expressing ###

for (group_var in group_vars) {
  filename <- paste0("dot_plot_", tolower(group_var), ".pdf")
  filepath <- file.path(output_dir, filename)
  
  if (!file.exists(filepath)) {
    message("Generating and saving: ", filepath)
    
    p <- DotPlot(
      object = seurat_obj_brain,
      features = genes_of_interest,
      group.by = group_var,
      cols = c("grey", "blue")
    ) +
      labs(x = "Gene", y = "Cell Type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plot_data <- p$data  # get gene and group labels here
    
    annots <- merge(
      plot_data,
      sig_df[, .(Gene, Lineage, SigLabel)],
      by.x = c("features.plot", "id"),
      by.y = c("Gene", "Lineage"),
      all.x = TRUE
    )
    
    annots$SigLabel[is.na(annots$SigLabel)] <- ""
    
    p <- p + geom_text(
      data = annots,
      aes(x = features.plot, y = id, label = SigLabel),
      size = 5,
      vjust = -0.6
    )
    
    
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

lineages <- c("All", "OPC", "ExNeu", "AST", "OL", "IN", "GLIALPROG")
age_levels <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years",
                "2-4 years", "4-10 years", "10-20 years", "Adult")
seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_levels)

# helper function to get summary data
get_summary_df <- function(seurat_obj, gene, lineage, age_levels) {
  if (lineage == "All") {
    seurat_sub <- seurat_obj
  } else {
    cells <- WhichCells(seurat_obj, expression = Lineage == lineage)
    seurat_sub <- subset(seurat_obj, cells = cells)
  }
  
  exprs_data <- FetchData(seurat_sub, vars = c(gene, "Age_Range"))
  
  summary_df <- exprs_data %>%
    mutate(Age_Range = factor(Age_Range, levels = age_levels)) %>%
    group_by(Age_Range) %>%
    summarize(
      avg_exprs = mean(.data[[gene]], na.rm = TRUE),
      percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100,
      .groups = "drop"
    ) %>%
    complete(Age_Range = age_levels, fill = list(avg_exprs = 0, percent_exprs = 0)) %>%
    mutate(Gene = gene, Lineage = lineage)
  
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

# find top 3 max avg. expressions
top_expr_rows <- combined_summary %>%
  arrange(desc(avg_exprs)) %>%
  slice(1:3)

excluded_gene <- top_expr_rows$Gene[1]
excluded_lineage <- top_expr_rows$Lineage[1]

# remove the top max row for consistent y-scale
target_summary <- combined_summary %>%
  filter(!(Gene == excluded_gene & Lineage == excluded_lineage))

# global axis scales
y_max <- ceiling(max(target_summary$avg_exprs, na.rm = TRUE) * 10) / 10
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
    
    df <- df %>% mutate(Age_Range = factor(Age_Range, levels = age_levels))
    
    p <- ggplot(df, aes(x = Age_Range, y = avg_exprs, group = 1)) +
      geom_line(color = "#2c7fb8", linewidth = 1) +
      geom_point(aes(size = percent_exprs), color = "#2c7fb8") +
      scale_size_continuous(limits = c(0, max_percent), breaks = scales::pretty_breaks(n = 5)) +
      theme_minimal() +
      labs(
        title = paste("Expression of", gene, "across Age Ranges"),
        subtitle = subtitle_text,
        x = "Age Range",
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
    if (gene == excluded_gene && lineage == excluded_lineage) {
      message("Skipping plot for max expression (", gene, ", ", lineage, ")")
      next
    }
    
    df <- target_summary %>%
      filter(Gene == gene, Lineage == lineage)
    
    generate_and_save_age_expr_plot(df, gene, lineage, y_max, max_percent, output_dir)
  }
}

genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")
celltype_map <- list(SORCS1 = "OPC", SORCS2 = "AST", SORCS3 = "OPC")

# Filter combined_summary to include only relevant gene-celltype combos
overlay_data <- combined_summary %>%
  filter(Gene %in% genes_of_interest) %>%
  rowwise() %>%
  filter(Lineage == celltype_map[[Gene]]) %>%
  ungroup() %>%
  select(Gene, Lineage, Age_Range, avg_exprs, percent_exprs)

# Save to CSV
overlay_path <- file.path(output_dir, "overlay_expression_summary.csv")
write.csv(overlay_data, overlay_path, row.names = FALSE)


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
age_counts <- as.data.frame(table(seurat_obj$Age_Range))
colnames(age_counts) <- c("Age_Range", "Cell_Count")

p2 <- ggplot(age_counts, aes(x = Age_Range, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Age Range", y = "Number of Cells", title = "Cell Counts per Age Range") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p2, "cell_counts_per_age_range.pdf")

### 4c. cell counts per broad region ###
region_counts <- as.data.frame(table(seurat_obj$Region_Broad))
colnames(region_counts) <- c("Region_Broad", "Cell_Count")

p3 <- ggplot(region_counts, aes(x = Region_Broad, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Broad Region", y = "Number of Cells", title = "Cell Counts per Broad Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p3, "cell_counts_per_broad_region.pdf")

### 4d. cell counts per lineage ###
lineage_counts <- as.data.frame(table(seurat_obj$Lineage))
colnames(lineage_counts) <- c("Lineage", "Cell_Count")

p4 <- ggplot(lineage_counts, aes(x = Lineage, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Lineage", y = "Number of Cells", title = "Cell Counts per Lineage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p4, "cell_counts_per_lineage.pdf")

### 4e. cell type composition per region (%) ###
region_celltype_counts <- seurat_obj@meta.data %>%
  group_by(Region_Broad, Lineage) %>%
  summarise(Count = n(), .groups = "drop")

region_celltype_percent <- region_celltype_counts %>%
  group_by(Region_Broad) %>%
  mutate(Percent = Count / sum(Count) * 100)

p5 <- ggplot(region_celltype_percent, aes(x = Region_Broad, y = Percent, fill = Lineage)) +
  geom_bar(stat = "identity") +
  labs(x = "Region", y = "Cell Type Composition (%)", fill = "Cell Type",
       title = "Cell Type Composition per Brain Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

generate_and_save_data_stats_plot(p5, "cell_type_composition_per_region.pdf", width = 10)


#### 5. Heatmap Analysis of New and Known Autism Risk Genes ####

# Load expression data and metadata
exprs_data_asd <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
meta <- seurat_obj@meta.data
meta$cell <- rownames(meta)

# Load rank data
rank_data <- readRDS(file.path(base_dir, "cell_rankings_velmeshev.rds"))

# Load ASD risk genes
setwd("C:/Abbie/research/seurat")
asd_risk_genes <- read.csv("./autism_risk_genes_combined.csv", sep = ",", header = TRUE)
asd_risk_genes <- asd_risk_genes[asd_risk_genes$Gene %in% rownames(seurat_obj), ]
gene_status <- setNames(asd_risk_genes$Status, asd_risk_genes$Gene)

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
    left_join(meta[, c("cell", "Region_Broad")], by = "cell")
  
  avg_exprs_region <- df_region %>%
    group_by(gene, Region_Broad) %>%
    summarise(avg_exprs = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = Region_Broad, values_from = avg_exprs)
  
  write_xlsx(avg_exprs_region, file_exprs_region)
} else {
  message("Average expression table by Region exists, loading...")
  avg_exprs_region <- readxl::read_xlsx(file_exprs_region)
}

### 5b. avg expression per gene by lineage ###
if (!file.exists(file_exprs_lineage)) {
  message("Generating average expression table by Lineage")
  
  df_lineage <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Lineage")], by = "cell")
  
  avg_exprs_lineage <- df_lineage %>%
    group_by(gene, Lineage) %>%
    summarise(avg_exprs = mean(expression), .groups = "drop") %>%
    pivot_wider(names_from = Lineage, values_from = avg_exprs)
  
  write_xlsx(avg_exprs_lineage, file_exprs_lineage)
} else {
  message("Average expression table by Lineage exists, loading...")
  avg_exprs_lineage <- readxl::read_xlsx(file_exprs_lineage)
}

### 5c. % expression per gene by region ###
if (!file.exists(file_pct_exprs_region)) {
  message("Generating % expression table by Region")
  
  df_region <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Region_Broad")], by = "cell")
  
  pct_exprs_region <- df_region %>%
    group_by(gene, Region_Broad) %>%
    summarise(pct_exprs = sum(expression > 0) / n() * 100, .groups = "drop") %>%
    pivot_wider(names_from = Region_Broad, values_from = pct_exprs)
  
  write_xlsx(pct_exprs_region, file_pct_exprs_region)
} else {
  message("% Expression table by Region exists, loading...")
  pct_exprs_region <- readxl::read_xlsx(file_pct_exprs_region)
}

### 5d. % expression per gene by lineage ###
if (!file.exists(file_pct_exprs_lineage)) {
  message("Generating % expression table by Lineage")
  
  df_lineage <- exprs_data_asd %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta[, c("cell", "Lineage")], by = "cell")
  
  pct_exprs_lineage <- df_lineage %>%
    group_by(gene, Lineage) %>%
    summarise(pct_exprs = sum(expression > 0) / n() * 100, .groups = "drop") %>%
    pivot_wider(names_from = Lineage, values_from = pct_exprs)
  
  write_xlsx(pct_exprs_lineage, file_pct_exprs_lineage)
} else {
  message("% Expression table by Lineage exists, loading...")
  pct_exprs_lineage <- readxl::read_xlsx(file_pct_exprs_lineage)
}

### 5e. avg rank per gene by region ###

if (!file.exists(file_rank_region)) {
  message("Generating average rank table by Region")
  
  df_region <- rank_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "rank") %>%
    left_join(meta[, c("cell", "Region_Broad")], by = "cell")
  
  avg_rank_region <- df_region %>%
    group_by(gene, Region_Broad) %>%
    summarise(avg_rank = mean(rank), .groups = "drop") %>%
    pivot_wider(names_from = Region_Broad, values_from = avg_rank)
  
  write_xlsx(avg_rank_region, file_rank_region)
} else {
  message("Average rank table by Region exists, loading...")
  avg_rank_region <- readxl::read_xlsx(file_rank_region)
}

### 5f. avg rank per gene by lineage ###

if (!file.exists(file_rank_lineage)) {
  message("Generating average rank table by Lineage")
  
  df_lineage <- rank_data %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "rank") %>%
    left_join(meta[, c("cell", "Lineage")], by = "cell")
  
  avg_rank_lineage <- df_lineage %>%
    group_by(gene, Lineage) %>%
    summarise(avg_rank = mean(rank), .groups = "drop") %>%
    pivot_wider(names_from = Lineage, values_from = avg_rank)
  
  write_xlsx(avg_rank_lineage, file_rank_lineage)
} else {
  message("Average rank table by Lineage exists, loading...")
  avg_rank_lineage <- readxl::read_xlsx(file_rank_lineage)
}

### 5g. heatmap by lineage ###
if (!file.exists(file_heatmap_lineage)) {
  message("Generating heatmap by Lineage")
  
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
    column_names_gp = gpar(fontsize = 4),
    col = exprs_colors,
    top_annotation = col_annot
  )
  
  pdf(file_heatmap_lineage, width = 20, height = 8)
  draw(ht_lineage)
  dev.off()
} else {
  message("Lineage heatmap file exists, skipping generation.")
}

### 5h. heatmap by region ###
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

### 5i. heatmap of % exprs by lineage ###
if (!file.exists(file_heatmap_pct_lineage)) {
  message("Generating % exprs heatmap by Lineage")
  
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
  message("% Exprs Lineage heatmap exists, skipping generation.")
}

### 5j. heatmap of % exprs by region ###
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


### 5k. heatmap for rank by region ###
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

### 5l. heatmap for rank by lineage ###
if (!file.exists(file_heatmap_rank_lineage)) {
  message("Generating rank heatmap by Lineage")
  
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
  message("Lineage rank heatmap file exists, skipping generation.")
}

