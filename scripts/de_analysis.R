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
library(dplyr)

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_common_genes.rds")
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj$celltype_age <- paste(seurat_obj$Lineage, seurat_obj$Age_Range, sep = "_")

# Define your Seurat object
Idents(seurat_obj) <- "celltype_age"

# Set your genes of interest and matching cell types
gene_celltype_map <- list(
  "SORCS1" = "OPC",
  "SORCS3" = "OPC",
  "SORCS2" = "AST"
)

# Extract the age ranges in order
all_groups <- unique(seurat_obj$celltype_age)
all_ages <- unique(gsub(".*_", "", all_groups))
# Manually define the age order (from your data)
age_order <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", 
               "4-10 years", "10-20 years", "Adult")


#### 1. DE analysis ####
# Create output list to store results
de_results <- list()

# Loop over each gene and its corresponding cell type
for (gene in names(gene_celltype_map)) {
  cell_type <- gene_celltype_map[[gene]]
  
  # Filter age stages that have that cell type
  age_groups <- paste0(cell_type, "_", age_order)
  existing <- age_groups[age_groups %in% all_groups]
  
  # Loop over adjacent age ranges
  for (i in 1:(length(existing) - 1)) {
    ident.1 <- existing[i + 1]
    ident.2 <- existing[i]
    
    message("Comparing ", ident.1, " vs ", ident.2, " for gene ", gene)
    
    # Try block in case some comparisons fail (e.g., too few cells)
    try({
      res <- FindMarkers(
        seurat_obj,
        ident.1 = ident.1,
        ident.2 = ident.2,
        features = gene,
        test.use = "wilcox",
        logfc.threshold = 0,
        verbose = FALSE
      )
      comp_name <- paste(gene, ident.1, "vs", ident.2, sep = "_")
      de_results[[comp_name]] <- res
    }, silent = TRUE)
  }
}

# View results
summary_list <- list()

# Loop over each result
for (name in names(de_results)) {
  df <- de_results[[name]]
  
  # Only proceed if the result is not empty
  if (!is.null(df) && nrow(df) > 0) {
    df$gene <- rownames(df)        # Add gene name column
    df$comparison <- name          # Add comparison name column
    rownames(df) <- NULL           # Reset row names
    summary_list[[length(summary_list) + 1]] <- df
  }
}

# Combine all into one data frame
summary_df <- do.call(rbind, summary_list)

# Save result
write.csv(summary_df, file.path(base_dir, "age_DE_summary.csv"), row.names = FALSE)


#### 2. Violin plot ####
for (gene in names(gene_celltype_map)) {
  cell_type <- gene_celltype_map[[gene]]
  
  # Subset object to relevant cell type
  subset_obj <- subset(seurat_obj, subset = Lineage == cell_type)
  
  # Ensure Age_Range is a factor in the right order
  subset_obj$Age_Range <- factor(subset_obj$Age_Range, levels = age_order)
  
  # Plot
  p <- VlnPlot(
    subset_obj,
    features = gene,
    group.by = "Age_Range",
    pt.size = 0
  ) + ggtitle(paste(gene, "expression in", cell_type, "across ages")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save
  pdf_path <- file.path(output_dir, paste0("violin_", gene, "_", cell_type, ".pdf"))
  ggsave(pdf_path, plot = p, width = 8, height = 5)
}

#### 3. Overlay line plot on violin plot ####

overlay_path <- file.path(output_dir, "overlay_expression_summary.csv")
overlay_data <- read.csv(overlay_path)

seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_order)
overlay_data$Age_Range <- factor(overlay_data$Age_Range, levels = age_order)


# Loop through genes and plot
for (gene in names(gene_celltype_map)) {
  cell_type <- gene_celltype_map[[gene]]
  
  # Subset Seurat object
  subset_obj <- subset(seurat_obj, subset = Lineage == cell_type)
  valid_cells <- which(!is.na(FetchData(subset_obj, vars = gene)[,1]))
  subset_obj <- subset_obj[, valid_cells]
  
  # Make sure Age_Range is a factor with proper order
  subset_obj$Age_Range <- factor(subset_obj$Age_Range, levels = age_order)
  
  # Subset overlay data
  df <- overlay_data %>%
    filter(Gene == gene, Lineage == cell_type)
  
  # Determine y-axis limits
  y_max <- max(c(df$avg_exprs, FetchData(subset_obj, vars = gene)[,1]), na.rm = TRUE) * 1.1
  max_percent <- max(df$percent_exprs, na.rm = TRUE)
  
  # Violin plot
  vln <- VlnPlot(
    subset_obj,
    features = gene,
    group.by = "Age_Range",
    pt.size = 0
  ) +
    ggtitle(paste(gene, "expression in", cell_type, "across ages")) +
    labs(x = "Age Range")
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Overlay line + dots
  overlay_plot <- vln +
    geom_line(
      data = df,
      aes(x = Age_Range, y = avg_exprs, group = 1),
      color = "black",
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = df,
      aes(x = Age_Range, y = avg_exprs, size = percent_exprs),
      color = "black",
      inherit.aes = FALSE
    ) +
    scale_size_continuous(name = "% Expressing", limits = c(0, max_percent)) + 
    coord_cartesian(ylim = c(0, y_max))

  
  # Save plot
  filename <- paste0("overlay_violin_", gene, "_", cell_type, ".pdf")
  ggsave(
    filename = file.path(output_dir, filename),
    plot = overlay_plot,
    width = 8,
    height = 10
  )
  
  message("Saved: ", filename)
}


