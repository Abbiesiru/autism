base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/prepostnatal",
                   "gauss" = "/home/abbiew/single_cell/velmeshev",
                   "macbook-air.lan"     = "/Users/abbiesiru/Desktop/research/shen lab/autism",
                   "."
)

output_folder <- "output/velmeshev/age_specificity"
output_dir <- file.path(base_dir, output_folder)

# Create folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(Seurat)
library(dplyr)
library(tibble)
library(readxl)

seurat_obj_path <- file.path(base_dir, "files/velmeshev/seurat_obj_subset_common_genes.rds")
seurat_obj <- readRDS(seurat_obj_path)
asd_risk_genes <- read.csv(file.path(base_dir, "files/autism_risk_genes_combined.csv"), sep = ",", header = TRUE)
asd_risk_genes <- asd_risk_genes[asd_risk_genes$Gene %in% rownames(seurat_obj), ]
asd_risk_genes <- asd_risk_genes$Gene
age_ranges <- unique(seurat_obj$Age_Range)

#### 1. generate matrices/marker gene peak list ####
binary_matrix <- matrix(0, nrow = length(asd_risk_genes), ncol = length(age_ranges))
rownames(binary_matrix) <- asd_risk_genes
colnames(binary_matrix) <- age_ranges


logp_matrix <- matrix(0, nrow = length(asd_risk_genes), ncol = length(celltypes),
                      dimnames = list(asd_risk_genes, celltypes))

gene_peak_list <- list()

for (age in age_ranges) {
  message("Processing: ", age)
  
  markers <- FindMarkers(
    seurat_obj,
    ident.1 = age,
    group.by = "Age_Range",
    only.pos = FALSE,
    logfc.threshold = 0,
    min.pct = 0.1,
    test.use = "wilcox"
  )
  
  # filter ASD genes and apply significance threshold
  markers_asd <- markers %>%
    rownames_to_column("gene") %>%
    filter(gene %in% asd_risk_genes) %>%
    mutate(is_specific = ifelse(avg_log2FC > 1 & p_val_adj < 0.05, 1, 0)) %>%
    column_to_rownames("gene")
  
  gene_peak_list[[age]] <- markers_asd
}

# determine unique age with max expression if multiple significant
gene_peak_df <- data.frame()

# iterate through each ASD gene
for (gene in asd_risk_genes) {
  sig_ages <- c()
  avg_expr <- c()
  logfc <- c()
  
  # loop through each age range
  for (age in age_ranges) {
    markers_age <- gene_peak_list[[age]]
    
    if (!is.null(markers_age) && gene %in% rownames(markers_age)) {
      m <- markers_age[gene, , drop = FALSE]
      
      if (!is.na(m$is_specific) && m$is_specific == 1) {
        sig_ages <- c(sig_ages, age)
        expr_vals <- FetchData(seurat_obj, vars = gene)[seurat_obj$Age_Range == age, 1]
        avg_expr <- c(avg_expr, mean(expr_vals, na.rm = TRUE))
        logfc <- c(logfc, m$avg_log2FC)
      }
    }
  }
  
  # if gene is significant in ≥1 age range
  if (length(sig_ages) > 0) {
    peak_idx <- which.max(avg_expr)
    peak_age <- sig_ages[peak_idx]
    
    # double-check dimensions
    if (gene %in% rownames(binary_matrix) && peak_age %in% colnames(binary_matrix)) {
      binary_matrix[gene, peak_age] <- 1
    }
    
    gene_peak_df <- rbind(
      gene_peak_df,
      data.frame(
        Gene = gene,
        Peak_Age = peak_age,
        Mean_Expr = avg_expr[peak_idx],
        Avg_Log2FC = logfc[peak_idx],
        stringsAsFactors = FALSE
      )
    )
    
    message("✅ ", gene, " → ", peak_age)
  }
}

# save outputs

saveRDS(gene_peak_list, file = file.path(output_dir, "gene_peak_list.rds"))


write.csv(as.data.frame(binary_matrix),
          file.path(output_dir, "ASD_gene_specificity_by_age_matrix.csv"))

write.csv(gene_peak_df,
          file.path(output_dir, "ASD_gene_peak_age.csv"),
          row.names = FALSE)


#### 2. heatmap ####
library(ComplexHeatmap)
library(circlize)
library(grid)

gene_peak_list <- readRDS(file.path(output_dir, "gene_peak_list.rds"))

# Create logp matrix
logp_matrix <- matrix(0, nrow = length(asd_risk_genes), ncol = length(age_ranges),
                      dimnames = list(asd_risk_genes, age_ranges))

for (age in age_ranges) {
  markers_asd <- gene_peak_list[[age]]
  if (!is.null(markers_asd)) {
    vals <- markers_asd$p_val_adj
    vals[vals == 0] <- 1e-300  # prevent -Inf
    logp_matrix[rownames(markers_asd), age] <- ifelse(markers_asd$is_specific == 1,
                                                      -log10(vals), 0)
  }
}


# order age ranges
age_order <- c(
  "2nd trimester",
  "3rd trimester",
  "0-1 years",
  "1-2 years",
  "2-4 years",
  "4-10 years",
  "10-20 years",
  "Adult"
)

age_order <- intersect(age_order, colnames(logp_matrix))  # keep only ages in your dataset
logp_matrix <- logp_matrix[, age_order, drop = FALSE]

cap_val <- quantile(logp_matrix[logp_matrix > 0], 0.99, na.rm = TRUE) # cap values to prevent inf
logp_matrix[logp_matrix > cap_val] <- cap_val

# color scale
grey_threshold <- -log10(0.05)
col_fun <- colorRamp2(c(0, grey_threshold, cap_val / 2, cap_val),
                      colors = c("grey80", "cadetblue1", "mediumpurple3", "plum4"))

pdf_width <- max(8, ncol(logp_matrix) * 1.2)
pdf_height <- max(8, nrow(logp_matrix) / 10)
pdf(file.path(output_dir, "ASD_gene_age_specificity_heatmap.pdf"),
    width = pdf_width, height = pdf_height)

# heatmap
Heatmap(logp_matrix,
        name = "-log10(adj p)",
        col = col_fun,
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        row_dend_side = "left",
        column_dend_side = "top",
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 12)
)

dev.off()


#### 3. correlation between effect size and expression peak age range ####

library(dplyr)

gene_peak_df <- read.csv(file.path(output_dir, "ASD_gene_peak_age.csv"), 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

gene_effect_sizes <- read_xlsx(file.path(base_dir, "files/gene_effect_sizes.xlsx"))
gene_effect_sizes$LoF_fold <- as.numeric(gene_effect_sizes$LoF_fold)
asd_effects <- gene_effect_sizes %>%
  filter(Gene %in% asd_risk_genes)

# take max LoF_fold per gene
asd_effects_max <- asd_effects %>%
  group_by(Gene) %>%
  summarise(
    Max_LoF_fold = {
      vals <- LoF_fold[!is.na(LoF_fold) & !is.infinite(LoF_fold)]
      if(length(vals) == 0) NA_real_ else max(vals)
    }
  ) %>%
  ungroup()

gene_peak_sub <- gene_peak_df %>%
  filter(Gene %in% asd_effects_max$Gene)
gene_peak_lof <- merge(gene_peak_sub, asd_effects_max, by = "Gene") # merge


# Assign numeric codes to age ranges (optional)
age_order <- c(
  "2nd trimester",
  "3rd trimester",
  "0-1 years",
  "1-2 years",
  "2-4 years",
  "4-10 years",
  "10-20 years",
  "Adult"
)
gene_peak_lof$Peak_Age <- factor(gene_peak_lof$Peak_Age, levels = age_order)
gene_peak_lof$Peak_Age_Num <- as.numeric(gene_peak_lof$Peak_Age)

# Spearman correlation (for ordinal variable)
gene_peak_lof$Max_LoF_fold <- as.numeric(gene_peak_lof$Max_LoF_fold)
cor_test <- cor.test(
  gene_peak_lof$Peak_Age_Num[!is.na(gene_peak_lof$Max_LoF_fold)],
  gene_peak_lof$Max_LoF_fold[!is.na(gene_peak_lof$Max_LoF_fold)],
  method = "spearman"
)

library(ggplot2)

# Ensure Max_LoF_fold is numeric and Peak_Age is a factor
gene_peak_lof$Max_LoF_fold <- as.numeric(gene_peak_lof$Max_LoF_fold)
gene_peak_lof$Peak_Age <- factor(gene_peak_lof$Peak_Age,
                                 levels = c(
                                   "2nd trimester",
                                   "3rd trimester",
                                   "0-1 years",
                                   "1-2 years",
                                   "2-4 years",
                                   "4-10 years",
                                   "10-20 years",
                                   "Adult"
                                 )
                                 )

# Remove rows where Max_LoF_fold is NA
plot_data <- gene_peak_lof[!is.na(gene_peak_lof$Max_LoF_fold), ]

# Create boxplot
p <- ggplot(plot_data, aes(x = Peak_Age, y = Max_LoF_fold)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(x = "Peak Expression Age",
       y = "Max LoF_fold",
       title = "ASD Gene Effect Size by Peak Expression Age") +
  theme_minimal()

# Save as PDF
pdf(file.path(output_dir, "ASD_gene_effect_size_vs_peak_age_boxplot.pdf"),
    width = 8, height = 6)
print(p)
dev.off()


library(ggplot2)


# Create violin plot
p <- ggplot(plot_data, aes(x = Peak_Age, y = Max_LoF_fold, fill = Peak_Age)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  labs(
    title = "ASD Gene Effect Size (LoF_fold) by Peak Expression Age",
    x = "Peak Expression Age",
    y = "Max LoF_fold"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# Save as PDF
pdf(file.path(output_dir, "ASD_gene_LoF_vs_PeakAge_violin.pdf"),
    width = 8, height = 6)
print(p)
dev.off()

# Save as PNG
png(file.path(output_dir, "ASD_gene_effect_size_vs_peak_age_violin.png"),
    width = 1200, height = 900, res = 150)
print(p)
dev.off()


