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
library(dplyr)
library(msigdbr)
library(fgsea)
library(ggplot2)

seurat_obj_path <- file.path(base_dir, "seurat_obj.rds")
seurat_obj <- readRDS(seurat_obj_path)


#### 1. Subset data to astrocytes and postnatal cells ####

astrocytes_postnatal <- subset(seurat_obj, subset = Lineage == "AST" 
                               & (Age_Range == "0-1 years" | Age_Range == "1-2 years" | 
                                    Age_Range == "2-4 years" | Age_Range == "4-10 years" | 
                                    Age_Range == "10-20 years" | Age_Range == "Adult"))


#### 2. Calculate Spearman rank correlation of each gene to SORCS1 ####

# Extract SORCS1 expression from matrix
expr_matrix <- GetAssayData(astrocytes_postnatal, slot = "data")
SORCS1_expr <- expr_matrix["SORCS1", ]

# Compute Spearman correlation between SORCS1 and all other genes
correlations <- apply(expr_matrix, 1, function(x) cor(x, SORCS1_expr, method = "spearman"))

# Compute p-values & adj p-values
pvals <- apply(expr_matrix, 1, function(x) cor.test(x, SORCS1_expr, method = "spearman")$p.value)
padj <- p.adjust(pvals, method = "fdr")

# Combine into data frame
cor_df <- data.frame(
  gene = rownames(expr_matrix),
  spearman_corr = correlations,
  pval = pvals,
  padj = padj
)


#### 3. GSEA using continuous Spearman rank ####

# Prepare ranked list
ranked_genes <- cor_df$spearman_corr
names(ranked_genes) <- cor_df$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Use MSigDB gene sets (e.g., GO Biological Processes)
msigdb_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

# Create list of gene sets
gene_sets <- split(msigdb_go$gene_symbol, msigdb_go$gs_name)

# Run fgsea
fgsea_res <- fgseaMultilevel(
  pathways = gene_sets,
  stats = ranked_genes,
  minSize = 15,
  maxSize = 500
)

fgsea_res$leadingEdge <- sapply(fgsea_res$leadingEdge, paste, collapse = ",")

# Save results
write.csv(cor_df, file.path(output_dir, "SORCS1_astrocytes_spearman.csv"), row.names = FALSE)
write.csv(fgsea_res, file.path(output_dir, "SORCS1_gsea_results.csv"), row.names = FALSE)

