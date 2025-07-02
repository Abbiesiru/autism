base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/prepostnatal",
                   "gauss" = "/home/abbiew/single_cell/velmeshev",
                   "."
)

library(dplyr)
library(Seurat)
library(patchwork)
library(Matrix)

meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)

mat <- readMM("matrix.mtx.gz")
barcodes <- read.delim("barcodes.tsv.gz", header = FALSE)
features <- read.delim("features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
rownames(mat) <- make.unique(features$V1) # gene names
colnames(mat) <- barcodes$V1 # cell barcodes

seurat_obj <- CreateSeuratObject(counts = mat, project = "prepostnatal", meta.data=meta)

umap <- read.delim("UMAP.coords.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
colnames(umap) <- c("barcode", "UMAP_1", "UMAP_2")
rownames(umap) <- umap$barcode
umap$barcode <- NULL
seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = "RNA")

seurat_obj <- NormalizeData(seurat_obj)

saveRDS(seurat_obj, "seurat_obj.rds")

### subset seurat_obj ###

genes_of_interest <- c("SORCS1", "SORCS2", "SORCS3")
ASD_risk_genes <- read.csv("autism_risk_genes_combined.csv", sep = ",", header = TRUE)

seurat_obj_path <- file.path(base_dir, "seurat_obj.rds")
seurat_obj <- readRDS(seurat_obj_path)

seurat_obj_subset <- seurat_obj[ASD_risk_genes$Gene, ]

# sort age ranges from earliest to latest
age_levels <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")
seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_levels)

saveRDS(seurat_obj_subset, "seurat_obj_subset.rds")




#### 2. Identify & Visualize spaciotemporal expression pattern of genes ####

### Feature Plot ###

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
        ggtitle(paste("Age Range:", age, "-", group_var, "-", gene))
    })
    
    return(wrap_plots(plots, ncol = 4))
    
  } else {
    # Single UMAP, not split by age
    FeaturePlot(seurat_obj, reduction = "umap", features = gene) +
      geom_text(
        data = centroids,
        aes(x = x, y = y, label = Group),
        size = 3,
        color = "black",
        inherit.aes = FALSE
      ) +
      ggtitle(paste("All Ages -", group_var, "-", gene))
  }
}

# 8 plots, one per age range
generate_umap_plot(seurat_obj, group_var = "Lineage", gene = "SORCS3", split_by_age = TRUE)

# 1 combined plot
generate_umap_plot(seurat_obj, group_var = "Region_Broad", gene = "SORCS1", split_by_age = FALSE)


### Dot Plot ###

DotPlot(
  object = seurat_obj, 
  features = genes_list, 
  group.by = 'Region_Broad', 
  cols = c("red", "orange", "yellow", "green", "blue", "purple", "pink", "grey", "brown")
)
# FC-frontal/prefrontal cortex, CC-cingulate cortex, TC-temporal cortex, IC-insular cortex, MC-motor cortex, CTX-cortex

### Scatter Plot ###

gene <- "SORCS2"

exprs_data <- FetchData(seurat_obj, vars = c(gene, "Age_Range"))

summary_df <- exprs_data %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]]),
    percent_exprs = sum(.data[[gene]] > 0) / n() * 100
  )


ggplot(summary_df, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "by Age Range"),
    subtitle = "All Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Cells Expressed"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## OPC cells only ##

opc_cells <- WhichCells(seurat_obj, expression = Lineage == "OPC") # oligodendrocyte precursor cells
seurat_opc <- subset(seurat_obj, cells = opc_cells)

# Fetch expression and metadata for OPCs only
exprs_data_opc <- FetchData(seurat_opc, vars = c(gene, "Age_Range"))

# Summarize
summary_opc <- exprs_data_opc %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_opc, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "OPC Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## OUT cells only ##

out_cells <- WhichCells(seurat_obj, expression = Lineage == "OUT") # none of the other lineages
seurat_out <- subset(seurat_obj, cells = out_cells)

# Fetch expression and metadata for OUT cells only
exprs_data_out <- FetchData(seurat_out, vars = c(gene, "Age_Range"))

# Summarize
summary_out <- exprs_data_out %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_out, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "OUT Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ExNeu cells only ##

exneu_cells <- WhichCells(seurat_obj, expression = Lineage == "ExNeu") # excitatory neuron
seurat_exneu <- subset(seurat_obj, cells = exneu_cells)

# Fetch expression and metadata for ExNeu cells only
exprs_data_exneu <- FetchData(seurat_exneu, vars = c(gene, "Age_Range"))

# Summarize
summary_exneu <- exprs_data_exneu %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_exneu, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "ExNeu Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## AST cells only ##

ast_cells <- WhichCells(seurat_obj, expression = Lineage == "AST")
seurat_ast <- subset(seurat_obj, cells = ast_cells)

# Fetch expression and metadata for AST cells only
exprs_data_ast <- FetchData(seurat_ast, vars = c(gene, "Age_Range"))

# Summarize
summary_ast <- exprs_data_ast %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_ast, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "AST Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## OL cells only ##

ol_cells <- WhichCells(seurat_obj, expression = Lineage == "OL") # oligodendrocyte
seurat_ol <- subset(seurat_obj, cells = ol_cells)

# Fetch expression and metadata for OL cells only
exprs_data_ol <- FetchData(seurat_ol, vars = c(gene, "Age_Range"))

# Summarize
summary_ol <- exprs_data_ol %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_ol, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "OL Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## IN cells only ##

in_cells <- WhichCells(seurat_obj, expression = Lineage == "IN") # interneuron/inhibitory neuron
seurat_in <- subset(seurat_obj, cells = in_cells)

# Fetch expression and metadata for IN cells only
exprs_data_in <- FetchData(seurat_in, vars = c(gene, "Age_Range"))

# Summarize
summary_in <- exprs_data_in %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_in, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "IN Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## GLIALPROG cells only ##

glialprog_cells <- WhichCells(seurat_obj, expression = Lineage == "GLIALPROG") # glial progenitor
seurat_glialprog <- subset(seurat_obj, cells = glialprog_cells)

# Fetch expression and metadata for GLIALPROG cells only
exprs_data_glialprog <- FetchData(seurat_glialprog, vars = c(gene, "Age_Range"))

# Summarize
summary_glialprog <- exprs_data_glialprog %>%
  group_by(Age_Range) %>%
  summarize(
    avg_exprs = mean(.data[[gene]], na.rm = TRUE),
    percent_exprs = sum(.data[[gene]] > 0, na.rm = TRUE) / n() * 100
  )

# Plot
ggplot(summary_glialprog, aes(x = Age_Range, y = avg_exprs, size = percent_exprs)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Expression of", gene, "across Age Ranges"),
    subtitle = "GLIALPROG Cells",
    x = "Age Range",
    y = "Average Expression",
    size = "% Expressing Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


### dataset statistics ###


## cell counts per donor ##
donor_counts <- as.data.frame(table(seurat_obj$Individual))
colnames(donor_counts) <- c("Individual", "Cell_Count")

ggplot(donor_counts, aes(x = Individual, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Individual", y = "Number of Cells", title = "Cell Counts per Donor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## cell counts per age range ##
age_counts <- as.data.frame(table(seurat_obj$Age_Range))
colnames(age_counts) <- c("Age_Range", "Cell_Count")

ggplot(age_counts, aes(x = Age_Range, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Age Range", y = "Number of Cells", title = "Cell Counts per Age Range") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## cell counts per broad region ##
region_counts <- as.data.frame(table(seurat_obj$Region_Broad))
colnames(region_counts) <- c("Region_Broad", "Cell_Count")

ggplot(region_counts, aes(x = Region_Broad, y = Cell_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Broad Region", y = "Number of Cells", title = "Cell Counts per Broad Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## cell types in each region ##

region_celltype_counts <- seurat_obj@meta.data %>%
  group_by(Region_Broad, Lineage) %>%
  summarise(Count = n(), .groups = "drop")

region_celltype_percent <- region_celltype_counts %>%
  group_by(Region_Broad) %>%
  mutate(Percent = Count / sum(Count) * 100)

ggplot(region_celltype_percent, aes(x = Region_Broad, y = Percent, fill = Lineage)) +
  geom_bar(stat = "identity") +
  labs(x = "Region", y = "Cell Type Composition (%)", fill = "Cell Type",
       title = "Cell Type Composition per Brain Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### 3. Analysis of all 100 New Autism Risk Genes ####

exprs_data_all <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
exprs_data_asd <- exprs_data_all[autism_risk_genes_new$Gene[autism_risk_genes_new$Gene %in% rownames(exprs_data_all)], ]

# Combine with metadata
meta <- seurat_obj@meta.data
meta$cell <- rownames(meta)

### average expression per gene by region & lineage ###

# Make df: gene | cell | expression | region
df_region <- exprs_data_asd %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
  left_join(meta[, c("cell", "Region_Broad")], by = "cell")

# Compute average expression per gene per region
avg_exprs_region <- df_region %>%
  group_by(gene, Region_Broad) %>%
  summarise(avg_exprs = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = Region_Broad, values_from = avg_exprs)

# Make df: gene | cell | expression | cell type
df_lineage <- exprs_data_asd %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
  left_join(meta[, c("cell", "Lineage")], by = "cell")

# Compute average expression per gene per cell type
avg_exprs_lineage <- df_lineage %>%
  group_by(gene, Lineage) %>%
  summarise(avg_exprs = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = Lineage, values_from = avg_exprs)

# write_xlsx(avg_exprs_region, "avg_exprs_region.xlsx")
# write_xlsx(avg_exprs_lineage, "avg_exprs_lineage.xlsx")

### heatmaps and clustering by lineage ###

exprs_mat_lineage <- as.matrix(avg_exprs_lineage[ , -1])         
rownames(exprs_mat_lineage) <- avg_exprs_lineage$gene

# assign colors
exprs_colors <- colorRamp2(
  c(min(exprs_mat_lineage), median(exprs_mat_lineage), max(exprs_mat_lineage)),
  c("blue", "white", "red")
)

# column annotations
lineage_anno <- HeatmapAnnotation(
  Lineage = colnames(exprs_mat_lineage),
  col = list(Lineage = structure(
    rainbow(length(unique(colnames(exprs_mat_lineage)))),
    names = colnames(exprs_mat_lineage)
  ))
)

Heatmap(
  exprs_mat_lineage,
  name = "Average Expression Across Cell Types",
  column_title = "Cell Type",
  row_title = "Gene",
  top_annotation = lineage_anno,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  col = exprs_colors
)


### heatmaps and clustering by region ###

exprs_mat_region <- as.matrix(avg_exprs_region[ , -1])         
rownames(exprs_mat_region) <- avg_exprs_region$gene

# assign colors
exprs_colors <- colorRamp2(
  c(min(exprs_mat_region), median(exprs_mat_region), max(exprs_mat_region)),
  c("blue", "white", "red")
)

# column annotations
region_anno <- HeatmapAnnotation(
  Region = colnames(exprs_mat_region),
  col = list(Region = structure(
    rainbow(length(unique(colnames(exprs_mat_region)))),
    names = colnames(exprs_mat_region)
  ))
)

Heatmap(
  t(exprs_mat_region),
  name = "Avg. Expr.",
  column_title = "ASD Risk Genes",
  row_title = "Cell_Type",
  cluster_rows = TRUE,
  cluster_columns = TRUE,  
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 6),
  col = exprs_colors
)



