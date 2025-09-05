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
library(ggplot2)

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_common_genes.rds")
seurat_obj <- readRDS(seurat_obj_path)
genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")


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
      # Remove zeros for median, handle empty case
      median_exprs = {
        nonzero <- .data[[gene]][.data[[gene]] > 0]
        if (length(nonzero) == 0) 0 else median(nonzero, na.rm = TRUE)
      },
      percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100,
      .groups = "drop"
    ) %>%
    complete(Age_Range = age_levels, fill = list(median_exprs = 0, percent_exprs = 0)) %>%
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

celltype_map <- list(SORCS1 = "OPC", SORCS2 = "AST", SORCS3 = "OPC")

# Filter combined_summary to include only relevant gene-celltype combos
overlay_data <- combined_summary %>%
  filter(Gene %in% genes_of_interest) %>%
  rowwise() %>%
  filter(Lineage == celltype_map[[Gene]]) %>%
  ungroup() %>%
  select(Gene, Lineage, Age_Range, median_exprs, percent_exprs)

# Save to CSV
overlay_path <- file.path(output_dir, "overlay_median_expression_summary.csv")
write.csv(overlay_data, overlay_path, row.names = FALSE)
