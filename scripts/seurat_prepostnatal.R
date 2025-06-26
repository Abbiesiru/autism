setwd("C:/Abbie/research/seurat/prepostnatal")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(Matrix)
library(data.table)


#### 1. Setup the Seurat Object ####

# meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
# 
# mat <- readMM("matrix.mtx.gz")
# barcodes <- read.delim("barcodes.tsv.gz", header = FALSE)
# features <- read.delim("features.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
# rownames(mat) <- make.unique(features$V1) # gene names
# colnames(mat) <- barcodes$V1 # cell barcodes
# 
# seurat_obj <- CreateSeuratObject(counts = mat, project = "prepostnatal", meta.data=meta)
# 
# umap <- read.delim("UMAP.coords.tsv.gz", header = FALSE, stringsAsFactors = FALSE)
# colnames(umap) <- c("barcode", "UMAP_1", "UMAP_2")
# rownames(umap) <- umap$barcode
# umap$barcode <- NULL
# seurat_obj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap), key = "UMAP_", assay = "RNA")
# 
# 
# seurat_obj <- NormalizeData(seurat_obj)
# 
# saveRDS(seurat_obj, "seurat_obj.rds")

seurat_obj <- readRDS("seurat_obj.rds")

# DimPlot(seurat_obj, reduction = "umap", group.by = "Lineage")


#### 2. Identify & Visualize spaciotemporal expression pattern of genes ####

genes_list <- c("SORCS1", "SORCS2", "SORCS3")

age_levels <- c("2nd trimester", "3rd trimester", "0-1 years", "1-2 years", "2-4 years", "4-10 years", "10-20 years", "Adult")
seurat_obj$Age_Range <- factor(seurat_obj$Age_Range, levels = age_levels)

### Feature Plot ###

umap_df <- data.frame(
  x = seurat_obj@reductions$umap@cell.embeddings[, 1],
  y = seurat_obj@reductions$umap@cell.embeddings[, 2],
  Lineage = seurat_obj$Lineage
)

centroids <- umap_df %>%
  group_by(Lineage) %>%
  summarize(x = mean(x), y = mean(y), .groups = "drop")

plots <- lapply(unique(seurat_obj$Age_Range), function(age) {
  subset_obj <- subset(seurat_obj, subset = Age_Range == age)
  
  FeaturePlot(subset_obj, reduction = "umap", features = "SORCS1") +
    geom_text(
      data = centroids,
      aes(x = x, y = y, label = Lineage),
      size = 3,
      color = "black"
    ) +
    ggtitle(paste("Age Range:", age))
})

wrap_plots(plots, ncol = 2)


### Dot Plot ###

DotPlot(
  object = seurat_obj, 
  features = genes_list, 
  group.by = 'Age_Range', 
  cols = c("red", "orange", "yellow", "green", "blue", "purple", "pink", "grey", "brown")
  )

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

out_cells <- WhichCells(seurat_obj, expression = Lineage == "OUT") # 
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

exneu_cells <- WhichCells(seurat_obj, expression = Lineage == "ExNeu")
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

ol_cells <- WhichCells(seurat_obj, expression = Lineage == "OL")
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

in_cells <- WhichCells(seurat_obj, expression = Lineage == "IN")
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

glialprog_cells <- WhichCells(seurat_obj, expression = Lineage == "GLIALPROG")
seurat_glialprog <- subset(seurat_obj, cells = glialprog_cells)

# Fetch expression and metadata for AST cells only
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


