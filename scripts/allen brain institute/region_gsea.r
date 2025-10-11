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

all_genes <- 13121

### 1a. asd risk genes ###
tmp <- read.csv(file.path(base_dir, "files/autism_risk_genes_combined.csv"), sep = ",", header = TRUE)
asd_all <- tmp$Gene
asd_new <- tmp$Gene[tmp$Status == "New"]
asd_known <- tmp$Gene[tmp$Status == "Known"]

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

# ALL GENES - loop for all regions & ASD genes
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
  K <- length(region_genes)                          # number of upreg genes in region
  n <- length(asd_genes)                             # number of ASD genes
  overlap_genes <- intersect(region_genes, asd_genes) # actual overlapping genes
  k <- length(overlap_genes)                         # overlap size
  N <- universe_size                                 # total genes
  
  # hypergeometric test: P(X >= k)
  pval <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
  
  # fold enrichment
  expected <- (K / N) * n
  fold <- ifelse(expected > 0, k / expected, NA)
  
  return(data.frame(
    K = K, n = n, k = k, N = N,
    pval = pval,
    fold_enrichment = fold,
    overlap_genes = paste(overlap_genes, collapse = ";")
  ))
}

results_all <- do.call(rbind, lapply(names(region_upreg), function(region) {
  df <- enrichment_test(region_upreg[[region]], asd_all, length(all_genes))
  df$Region <- region
  return(df)
}))

# adjust p-values for multiple testing
### eigen value decomposition ###

all_genes <- unique(unlist(region_upreg))
region_names <- names(region_upreg)

region_gene_matrix <- matrix(0, nrow=length(region_names), ncol=length(all_genes),
                             dimnames=list(region_names, all_genes))

for (r in region_names) {
  region_gene_matrix[r, region_upreg[[r]]] <- 1
}

region_corr <- cor(t(region_gene_matrix))  # correlation between regions

eig <- eigen(region_corr, only.values = TRUE)$values
total_variance <- sum(eig)
proportion_variance_explained <- eig / total_variance

M <- nrow(region_corr)  # number of regions
# Meff <- which(cumsum(proportion_variance_explained) >= 0.95)[1]
# Meff <- sum(proportion_variance_explained>=0.01)

raw_pvals <- results_all$pval

# Bonferroni using Meff
bonf_corrected <- pmin(raw_pvals * Meff, 1)  # cap at 1
results_all$pval_adj <- bonf_corrected

# reorder columns
results_all <- results_all[, c("Region", "K", "n", "k", "overlap_genes", "N", "fold_enrichment", "pval", "FDR")]


#### 3. filter by significant results only ####
sig_results_all <- subset(results_all, FDR < 0.05)

# order by significance
sig_results_all <- sig_results_all[order(sig_results_all$FDR), ]









# NEW GENES - loop for all regions & ASD genes
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
  overlap_genes <- intersect(region_genes, asd_genes) # actual overlapping genes
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
    fold_enrichment = fold,
    overlap_genes = paste(overlap_genes, collapse = ";")
  ))
}

results_new <- do.call(rbind, lapply(names(region_upreg), function(region) {
  df <- enrichment_test(region_upreg[[region]], asd_new, length(all_genes))
  df$Region <- region
  return(df)
}))

# adjust p-values for multiple testing
raw_pvals <- results_new$pval
bonf_corrected <- pmin(raw_pvals * Meff, 1)  # cap at 1
results_new$FDR<- bonf_corrected

# reorder columns
results_new <- results_new[, c("Region", "K", "n", "k", "overlap_genes", "N", "fold_enrichment", "pval", "FDR")]


#### 3. filter by significant results only ####
sig_results_new <- subset(results_new, FDR < 0.05)

# order by significance
sig_results_new <- sig_results_new[order(sig_results_new$FDR), ]









# KNOWN GENES - loop for all regions & ASD genes
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
  overlap_genes <- intersect(region_genes, asd_genes) # actual overlapping genes
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
    fold_enrichment = fold,
    overlap_genes = paste(overlap_genes, collapse = ";")
  ))
}

results_known <- do.call(rbind, lapply(names(region_upreg), function(region) {
  df <- enrichment_test(region_upreg[[region]], asd_known, length(all_genes))
  df$Region <- region
  return(df)
}))

# adjust p-values for multiple testing
raw_pvals <- results_known$pval
bonf_corrected <- pmin(raw_pvals * Meff, 1)  # cap at 1
results_known$FDR<- bonf_corrected

# reorder columns
results_known <- results_known[, c("Region", "K", "n", "k", "overlap_genes", "N", "fold_enrichment", "pval", "FDR")]


#### 3. filter by significant results only ####
sig_results_known <- subset(results_known, FDR < 0.05)

# order by significance
sig_results_known <- sig_results_known[order(sig_results_known$FDR), ]











##### 4. MERGE

merged_results <- sig_results_all %>%
  # attach overlap_genes from the "new" subset
  left_join(
    results_new %>% 
      select(Region, overlap_genes) %>%
      rename(overlap_genes_new = overlap_genes),
    by = "Region"
  ) %>%
  # attach overlap_genes from the "known" subset
  left_join(
    results_known %>%
      select(Region, overlap_genes) %>%
      rename(overlap_genes_known = overlap_genes),
    by = "Region"
  ) %>%
  # add counts for new/known genes
  mutate(
    number_new = ifelse(is.na(overlap_genes_new) | overlap_genes_new == "",
                   0,
                   lengths(strsplit(overlap_genes_new, ";"))),
    number_known = ifelse(is.na(overlap_genes_known) | overlap_genes_known == "",
                     0,
                     lengths(strsplit(overlap_genes_known, ";")))
  ) %>%
  # rename and drop columns
  rename(
    region_size = K,
    overlap = k,
    universe_size = N
  ) %>%
  select(-n)   # drop column n

## write out
write.csv(results_all, file.path(output_dir, "asd_genes_by_region_all.csv"))
write.csv(results_new, file.path(output_dir, "asd_genes_by_region_new.csv"))
write.csv(results_known, file.path(output_dir, "asd_genes_by_region_known.csv"))
write.csv(merged_results, file.path(output_dir, "asd_genes_by_region_merged.csv"))

#### 5. DIRECT COMPARISON: known vs new genes per region ####

compare_known_new <- lapply(names(region_upreg), function(region) {
  region_genes <- region_upreg[[region]]
  
  # counts
  a <- sum(asd_known %in% region_genes)   # known in region
  b <- length(asd_known) - a              # known not in region
  c <- sum(asd_new %in% region_genes)     # new in region
  d <- length(asd_new) - c                # new not in region
  
  tab <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE,
                dimnames = list(GeneSet = c("Known","New"),
                                RegionUp = c("Yes","No")))
  
  ft <- fisher.test(tab)
  
  data.frame(
    Region = region,
    known_in_region = a,
    known_not_in_region = b,
    new_in_region = c,
    new_not_in_region = d,
    odds_ratio = unname(ft$estimate),
    pvalue_fisher = ft$p.value,
    stringsAsFactors = FALSE
  )
})

results_compare <- do.call(rbind, compare_known_new)
results_compare$FDR <- p.adjust(results_compare$pvalue_fisher, method = "BH")

# order by significance
results_compare <- results_compare[order(results_compare$pvalue_fisher), ]


### eigen value decomposition ###

all_genes <- unique(unlist(region_upreg))
region_names <- names(region_upreg)

region_gene_matrix <- matrix(0, nrow=length(region_names), ncol=length(all_genes),
                             dimnames=list(region_names, all_genes))

for (r in region_names) {
  region_gene_matrix[r, region_upreg[[r]]] <- 1
}

region_corr <- cor(t(region_gene_matrix))  # correlation between regions

eig <- eigen(region_corr, only.values = TRUE)$values
total_variance <- sum(eig)
proportion_variance_explained <- eig / total_variance

# Li & Ji method for effective number of independent tests
M <- nrow(region_corr)  # number of regions
# Meff <- sum(eig>=1)  # effective # of tests
Meff <- sum(proportion_variance_explained>=0.01)

raw_pvals <- results_compare$pvalue_fisher  # your list of Fisher test p-values

# Bonferroni using Meff
bonf_corrected <- pmin(raw_pvals * Meff, 1)  # cap at 1
# Custom BH-FDR with Meff
# bh_corrected <- p.adjust(raw_pvals, method="BH") * (length(raw_pvals) / Meff)
# bh_corrected <- pmin(bh_corrected, 1)  # cap at 1

results_compare$FDR_corrected <- bonf_corrected

# save
write.csv(results_compare, file.path(output_dir, "asd_known_vs_new_by_region.csv"), row.names = FALSE)


### scree plot ###
# w/ cumulative line 
cum_var <- cumsum(variance_explained)
plot(variance_explained, type = "b", pch = 19,
     xlab = "Component",
     ylab = "% Variance Explained",
     main = "Scree Plot of Brain Region Correlation Eigenvalues",
     ylim = c(0, max(cum_var)))  # extend y-axis to fit cumulative

lines(1:length(cum_var), cum_var, type = "b", col = "blue", pch = 17)
legend("topright", legend = c("Individual", "Cumulative"),
       col = c("black", "blue"), pch = c(19, 17), bty = "n")

# w/o
variance_explained <- eig / sum(eig) * 100  

# Scree plot showing % variance explained
plot(variance_explained, type = "b", pch = 19,
     xlab = "Component",
     ylab = "% Variance Explained",
     main = "Scree Plot of Brain Region Correlation Eigenvalues")


