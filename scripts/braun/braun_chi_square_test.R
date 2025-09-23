base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/braun",
                   "gauss" = "/home/abbiew/single_cell/braun",
                   "macbook-air.lan"     = "/Users/abbiesiru/Desktop/research/shen lab/autism",
                   "."
)

output_folder <- "output/braun"
output_dir <- file.path(base_dir, output_folder)

# Create folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(Seurat)
library(data.table)
library(dplyr)

seurat_obj_path <- file.path(base_dir, "files/braun/seurat_obj_merged_layers_only_processed.rds")
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj$Subclass <- recode(
  seurat_obj$Subclass,
  "VSMC" = "Vascular smooth muscle cells",
  "Schwann" = "Schwann cells",
  "PREAC" = "Pre-astrocytes",
  "OPC" = "Oligodendrocyte precursors",
  "Neuronal IPC" = "Neuronal intermediate progenitors",
  "Neuron" = "Neurons",
  "Neuroblast" = "Neuroblasts",
  "Immune" = "Immune cells",
  "Glioblast" = "Glioblasts",
  "Fibroblast" = "Fibroblasts",
  "Erythrocyte" = "Erythrocytes",
  "Endothelial" = "Endothelial cells",
  "COPs (premyelinating)" = "Committed oligodendrocyte progenitors"
)


#### CELL TYPE ####
# Lineage cell types
genes <- c("SORCS1", "SORCS2", "SORCS3")
cell_types <- c("Radial glia", 
                "Pre-astrocytes", 
                "Oligodendrocyte precursors", 
                "Committed oligodendrocyte progenitors", 
                "Neurons", 
                "Neuronal intermediate progenitors", 
                "Neuroblasts", 
                "Glioblasts")

# Get expression data and metadata
expr_data <- FetchData(seurat_obj, vars = genes)
subclass <- seurat_obj$Subclass

# Create one big data.table
df <- data.table(Cell = colnames(seurat_obj),
                 Subclass = subclass,
                 expr_data)

# Melt the expression data to long format (gene, expression)
df_long <- melt(df, id.vars = c("Cell", "Subclass"),
                variable.name = "Gene", value.name = "Expression")

# Binary expression: 1 = expressed, 0 = not expressed
df_long[, Expressed := as.integer(Expression > 0)]
df_long[, Expression := NULL]  # drop raw expression to save memory

# Preallocate results
results <- list()

# Loop through each gene and cell type
for (gene in genes) {
  for (ct in cell_types) {
    
    subset_dt <- df_long[Gene == gene]
    
    # Logical masks
    in_type <- subset_dt$Subclass == ct
    out_type <- !in_type
    
    expressed <- subset_dt$Expressed == 1
    
    # Use numeric to avoid integer overflow
    a <- as.numeric(sum(expressed & in_type))
    b <- as.numeric(sum(!expressed & in_type))
    c <- as.numeric(sum(expressed & out_type))
    d <- as.numeric(sum(!expressed & out_type))
    
    contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    # Safeguards: skip invalid tables
    if (any(rowSums(contingency) == 0) || any(colSums(contingency) == 0)) {
      raw_p <- NA
    } else {
      expected <- sum(contingency[1,]) * sum(contingency[,1]) / sum(contingency)
      
      raw_p <- tryCatch({
        if (a < expected) 1 else chisq.test(contingency)$p.value
      }, error = function(e) {
        NA
      })
    }
    
    results[[length(results) + 1]] <- list(
      Gene = gene,
      Subclass = ct,
      A = a, B = b, C = c, D = d,
      RawP = raw_p
    )
  }
}


# Combine and adjust p-values
results_df <- rbindlist(results)
results_df[, AdjustedP := pmin(RawP * .N, 1)]
results_df[, Significant := AdjustedP < 0.05]
results_df[, SigLabel := fifelse(is.na(AdjustedP), "",
                                 fifelse(AdjustedP < 0.001, "***",
                                         fifelse(AdjustedP < 0.01, "**",
                                                 fifelse(AdjustedP < 0.05, "*", ""))))]


# Save or view
fwrite(results_df, file.path(base_dir, "optimized_chisq_results_subclass.csv"))
print(results_df[Significant == TRUE])

#### BRAIN REGION ####
brain_regions <- unique(seurat_obj$Region_Broad)
regions <- seurat_obj$Region_Broad

# Melt to long format
df_long <- melt(df, id.vars = c("Cell", "Region"),
                variable.name = "Gene", value.name = "Expression")

# Binary expression status
df_long[, Expressed := as.integer(Expression > 0)]
df_long[, Expression := NULL]

# Initialize result list
results <- list()

# Loop through genes and brain regions
for (gene in genes) {
  for (region in brain_regions) {
    
    subset_dt <- df_long[Gene == gene]
    
    # Logical groups
    in_region <- subset_dt$Region == region
    out_region <- !in_region
    
    expressed <- subset_dt$Expressed == 1
    not_expressed <- !expressed
    
    # Counts for 2x2 table
    a <- sum(expressed & in_region)
    b <- sum(not_expressed & in_region)
    c <- sum(expressed & out_region)
    d <- sum(not_expressed & out_region)
    
    contingency <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    
    expected <- sum(contingency[1,]) * sum(contingency[,1]) / sum(contingency)
    
    raw_p <- if (a < expected) 1 else chisq.test(contingency)$p.value
    
    results[[length(results) + 1]] <- list(
      Gene = gene,
      Region = region,
      A = a, B = b, C = c, D = d,
      RawP = raw_p
    )
  }
}

# Combine and adjust p-values
results_df <- rbindlist(results)
results_df[, AdjustedP := pmin(RawP * .N, 1)]
results_df[, Significant := AdjustedP < 0.05]

# Save and view
fwrite(results_df, "optimized_chisq_results_region.csv")
print(results_df[Significant == TRUE])
