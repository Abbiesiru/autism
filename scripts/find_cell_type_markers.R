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
library(tibble)

setwd("C:/Abbie/research/seurat")
asd_risk_genes <- read.csv("./autism_risk_genes_combined.csv", sep = ",", header = TRUE)

seurat_obj_path <- file.path(base_dir, "seurat_obj_subset_common_genes.rds")
seurat_obj <- readRDS(seurat_obj_path)

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

asd_risk_genes <- asd_risk_genes[asd_risk_genes$Gene %in% rownames(seurat_obj), ]
asd_risk_genes <- asd_risk_genes$Gene

# Get unique cell types, excluding Outliers
celltypes <- setdiff(unique(seurat_obj$Lineage), "Outliers")

# Initialize binary matrix
binary_matrix <- matrix(0, nrow = length(asd_risk_genes), ncol = length(celltypes))
rownames(binary_matrix) <- asd_risk_genes
colnames(binary_matrix) <- celltypes

# Loop through each cell type and compute markers
for (ct in celltypes) {
  message("Processing: ", ct)
  
  markers <- FindMarkers(
    seurat_obj,
    ident.1 = ct,
    group.by = "Lineage",
    only.pos = FALSE,
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  # Subset to ASD genes
  markers_asd <- markers %>%
    rownames_to_column("gene") %>%
    filter(gene %in% asd_risk_genes) %>%
    mutate(is_specific = ifelse(avg_log2FC > 1 & p_val_adj < 0.05, 1, 0)) %>%
    column_to_rownames("gene")
  
  # Fill in binary matrix
  binary_matrix[rownames(markers_asd), ct] <- markers_asd$is_specific
}

# Convert to data frame and write out
binary_df <- as.data.frame(binary_matrix)
write.csv(binary_df, file.path(output_dir, "ASD_gene_specificity_matrix.csv"))

#### 1. Identify cell type(s) with highest average expression for each ASD gene ####
gene_celltype_df <- data.frame(
  gene = rownames(binary_df),
  celltypes = apply(binary_df, 1, function(x) {
    paste(names(x)[x == 1], collapse = ", ")
  }),
  row.names = NULL
)
write.csv(gene_celltype_df, file.path(output_dir, "specificity_matrix_by_gene.csv"))


#### 2. Identify cell type with highest average expression for each ASD gene ####
# Only keep genes that are markers for at least 1 cell type
binary_df <- read.csv(file.path(output_dir, "ASD_gene_specificity_matrix.csv"), row.names = 1)
marker_genes <- rownames(binary_df)[rowSums(binary_df) > 0]

# Initialize matrix
max_expr_matrix <- matrix(0, nrow = length(marker_genes), ncol = length(celltypes))
rownames(max_expr_matrix) <- marker_genes
colnames(max_expr_matrix) <- celltypes

# Loop through genes
for (gene in marker_genes) {
  avg_exprs <- sapply(celltypes, function(ct) {
    cells <- colnames(seurat_obj)[seurat_obj$Lineage == ct]
    if (length(cells) == 0) return(NA)
    mean(FetchData(seurat_obj, vars = gene)[cells, 1], na.rm = TRUE)
  })
  
  # Identify cell type with max average expression
  max_ct <- names(avg_exprs)[which.max(avg_exprs)]
  
  # Fill in matrix
  max_expr_matrix[gene, max_ct] <- 1
}

# Convert to data frame and write out
max_expr_df <- as.data.frame(max_expr_matrix)
write.csv(max_expr_df, file.path(output_dir, "ASD_max_expression_matrix.csv"))


#### 3. Plot % of known vs novel ASD genes expressed in neurons ####
library(ggplot2)
library(dplyr)

# Make sure rownames are genes
binary_df$Gene <- rownames(binary_df)

# Merge with known/novel annotation
setwd("C:/Abbie/research/seurat")
gene_class_df <- read.csv("./autism_risk_genes_combined.csv", sep = ",", header = TRUE)
gene_class_df <- gene_class_df[, !names(gene_class_df) %in% "BF.Type"]
merged <- merge(binary_df, gene_class_df, by = "Gene")

### 2a. Ex/In lineage separated ###

# Define lineages
lineages <- c("Excitatory neurons", "Inhibitory neurons")

# Calculate percentage expressed in each lineage per group
results <- merged %>%
  select(all_of(lineages), Status) %>%
  group_by(Status) %>%
  summarise(across(everything(), ~ mean(. == 1) * 100)) %>%
  pivot_longer(cols = -Status, names_to = "lineage", values_to = "percent_expressed")

# Plot
p1 <- ggplot(results, aes(x = Status, y = percent_expressed, fill = Status)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  facet_wrap(~lineage, scales = "free_y") +
  labs(
    x = "ASD Gene Group",
    y = "Percent Expressed",
    fill = "Group",
    title = "Specific Expression of Known vs Novel ASD Genes in Neuronal Lineages"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 14))

filename <- "known_v_novel_ex_in_neurons.pdf"
ggsave(
  filename = file.path(output_dir, filename),
  plot = p1,
  width = 10,
  height = 10
)

### 2b. All neurons merged ###

merged$neuronal_lineage <- (merged$`Excitatory neurons` == 1 | merged$`Inhibitory neurons` == 1)

# Calculate percentage of each group expressed in neuronal lineage
plot_data <- merged %>%
  group_by(Status) %>%
  summarise(percent_expressed = mean(neuronal_lineage) * 100)

# Plot as bar chart
p2 <- ggplot(plot_data, aes(x = Status, y = percent_expressed, fill = Status)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(
    x = "ASD Gene Group",
    y = "Percent Expressed in Neuronal Lineage",
    title = "Specific Expression of Known vs Novel ASD Genes in Neuronal Lineage"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Known" = "#377eb8", "New" = "#e41a1c")) +
  theme(text = element_text(size = 14))

filename <- "known_v_novel_combined.pdf"
ggsave(
  filename = file.path(output_dir, filename),
  plot = p2,
  width = 10,
  height = 10
)

