base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/prepostnatal",
                   "gauss" = "/home/abbiew/single_cell/velmeshev",
                   "macbook-air.lan"     = "/Users/abbiesiru/Desktop/research/shen lab/autism",
                   "."
)

output_folder <- "output/region/gsea"
output_dir <- file.path(base_dir, output_folder)

# Create folder if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

library(Seurat)
library(dplyr)
library(tibble)
library(grid)  

#### 1. obtain gene lists ####

all_genes <- 20000

### 1a. asd risk genes ###
tmp <- read.csv(file.path(base_dir, "files/autism_risk_genes_combined.csv"), sep = ",", header = TRUE)
asd_risk_genes <- tmp$Gene

### 1b. upregulated genes by brain region ###
upreg_genes <- read.csv(
  file.path(base_dir, "files/allen brain institute/region_gene_upreg_list.txt"),
  header = FALSE,
  sep = "\t",
  row.names = NULL,
  check.names = FALSE
)

ncols <- ncol(upreg_genes)
region_upreg <- list()
current_reg <- NULL

# functions  for determining if a row is a continuation of the last row/region (blank), and 
# if the cell contains a gene or enrichment level (is numeric)
is_blank <- function(x) {
  if (is.null(x)) return(TRUE)
  if (is.na(x)) return(TRUE)
  if (!nzchar(trimws(as.character(x)))) return(TRUE)
  FALSE
}
is_numeric_only <- function(x) {
  grepl("^[-+]?[0-9]+(\\.[0-9]+)?$", trimws(as.character(x)))
}

# loop for all regions & genes
for (i in seq_len(nrow(upreg_genes))) {
  row <- upreg_genes[i, ]
  
  # determine if this is a new region row (col2 blank) or continuation (col2 filled)
  if (is_blank(row[[2]])) {
    # new region row: col1 is region name, genes start at col3
    current_reg <- as.character(row[[1]])
    if (is_blank(current_reg)) {
      warning(sprintf("Row %d: blank region name — skipping", i)); next
    }
    if (is.null(region_upreg[[current_reg]])) region_upreg[[current_reg]] <- character(0)
    start_col <- 3
  } else {
    # continuation row: genes start at col1
    if (is.null(current_reg)) stop("Found continuation row before any region header at row ", i)
    start_col <- 1
  }
  
  # process each cell individually 
  for (j in start_col:ncols) {
    cell <- as.character(row[[j]])
    if (is_blank(cell)) next
    
    # some cells might contain multiple entries separated by semicolon/pipe; split them first
    entries <- unlist(strsplit(cell, "[;|]"))
    for (entry in entries) {
      entry <- trimws(entry)
      if (entry == "") next
      
      # extract gene part: everything before first comma/colon/semicolon/pipe
      gene <- sub("[,;:|].*$", "", entry)
      gene <- trimws(gene)
      
      # skip if it's purely numeric (enrichment score split into its own cell)
      if (is_numeric_only(gene)) next
      if (gene == "") next
      
      # append (unique to avoid duplicates)
      region_upreg[[current_reg]] <- unique(c(region_upreg[[current_reg]], gene))
    }
  }
}

# check
lengths_by_region <- sapply(region_upreg, length)
print(head(lengths_by_region, 20))
str(region_upreg)  

#### 2. enrichment analysis ####
enrichment_test <- function(region_genes, asd_genes, universe_size) {
  K <- length(region_genes)                         # number of upreg genes in region
  n <- length(asd_genes)                            # number of ASD genes
  k <- sum(asd_genes %in% region_genes)             # overlap
  N <- universe_size                                # total genes
  
  # hypergeometric test: P(X >= k)
  pval <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
  
  # fold enrichment
  expected <- (K / N) * n
  fold <- ifelse(expected > 0, k / expected, NA)
  
  return(data.frame(
    K = K, n = n, k = k, N = N,
    pval = pval,
    fold_enrichment = fold
  ))
}

results <- do.call(rbind, lapply(names(region_upreg), function(region) {
  df <- enrichment_test(region_upreg[[region]], asd_risk_genes, all_genes)
  df$Region <- region
  return(df)
}))

# adjust p-values for multiple testing
results$FDR <- p.adjust(results$pval, method = "BH")

# reorder columns
results <- results[, c("Region", "K", "n", "k", "N", "fold_enrichment", "pval", "FDR")]


#### 3. filter by significant results only ####
sig_results <- subset(results, FDR < 0.05 & fold_enrichment > 1.5)

# order by significance
sig_results <- sig_results[order(sig_results$FDR), ]

