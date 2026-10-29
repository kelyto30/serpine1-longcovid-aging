# =============================================================================
# Script 01: Data Download, Quality Control, and Normalization
# Project:   Peripheral SERPINE1/PAI-1 as a Cross-Compartment Biomarker
#            of Long COVID-Associated Cognitive Aging
# Dataset:   GSE224615 (PBMC RNA-seq; Yin et al., 2023, Nature Immunology)
# Author:    Kelechi Wisdom Elechi
# PI:        Gabriel A. de Erausquin
# Affil.:    UT Health San Antonio / Barshop Institute
# Date:      2026-10
# =============================================================================
# INPUT:  GSE224615_DEGs.xlsx  (downloaded from GEO or placed in data/)
# OUTPUT: results/expr_mat.rds        -- log2(CPM+1) matrix (18,040 genes x 35 samples)
#         results/col_data.rds        -- sample metadata data.frame
#         results/figures/fig1_pca.png
#         results/figures/fig_library_size.png
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

# ── 1. Paths ──────────────────────────────────────────────────────────────────
DATA_DIR    <- "data"
RESULTS_DIR <- "results"
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,     showWarnings = FALSE, recursive = TRUE)

# GSE224615 count matrix (Excel format as distributed on GEO)
count_file <- file.path(DATA_DIR, "GSE224615_DEGs.xlsx")

if (!file.exists(count_file)) {
  stop(
    "Count file not found: ", count_file, "\n",
    "Download GSE224615_DEGs.xlsx from:\n",
    "  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224615\n",
    "and place it in the data/ directory."
  )
}

# ── 2. Load count matrix ──────────────────────────────────────────────────────
message("Loading count matrix from: ", count_file)
raw_xl <- read_excel(count_file)

# The first column is Gene.name (character); all remaining columns are samples
sample_cols <- setdiff(names(raw_xl), "Gene.name")
message("  Detected ", length(sample_cols), " sample columns")

# Extract gene names and numeric expression matrix
gene_names <- raw_xl[["Gene.name"]]
expr_mat_raw <- as.data.frame(raw_xl[, sample_cols])
expr_mat_raw <- apply(expr_mat_raw, 2, as.numeric)   # ensure numeric
rownames(expr_mat_raw) <- gene_names

# Drop any NA rows produced by the Gene.name column artefact
expr_mat_raw <- expr_mat_raw[!is.na(rownames(expr_mat_raw)), ]
message("  Expression matrix: ", nrow(expr_mat_raw), " genes x ",
        ncol(expr_mat_raw), " samples")

# ── 3. Build sample metadata ──────────────────────────────────────────────────
# Coding: sample IDs that map to Long COVID (LC) or non-LC status
# Based on Yin et al. 2023 supplementary metadata (GSE224615)
# LC samples: GSM7027459 - GSM7027480 (22 samples)
# NonLC samples: GSM7027481 - GSM7027493 (13 samples)
all_ids <- colnames(expr_mat_raw)

col_data <- data.frame(
  sample_id  = all_ids,
  row.names  = all_ids,
  stringsAsFactors = FALSE
)

# Assign group based on GSM range
gsm_num <- as.integer(sub("GSM", "", col_data$sample_id))
col_data$group <- ifelse(gsm_num <= 7027480, "LC", "NonLC")
col_data$group <- factor(col_data$group, levels = c("NonLC", "LC"))

message("  Sample groups:\n",
        paste(capture.output(table(col_data$group)), collapse = "\n"))

# ── 4. Normalization: log2(CPM + 1) ──────────────────────────────────────────
lib_sizes   <- colSums(expr_mat_raw, na.rm = TRUE)
cpm_mat     <- sweep(expr_mat_raw, 2, lib_sizes / 1e6, FUN = "/")
log2cpm_mat <- log2(cpm_mat + 1)

message("Normalization complete. Value range: [",
        round(min(log2cpm_mat, na.rm = TRUE), 2), ", ",
        round(max(log2cpm_mat, na.rm = TRUE), 2), "]")

# ── 5. QC: Library size distribution ─────────────────────────────────────────
lib_df <- data.frame(
  sample = names(lib_sizes),
  lib_size = lib_sizes,
  group = col_data[names(lib_sizes), "group"]
)

p_lib <- ggplot(lib_df, aes(x = reorder(sample, lib_size), y = lib_size / 1e6,
                              fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  coord_flip() +
  labs(
    title = "GSE224615: Library sizes per sample",
    x = NULL, y = "Total counts (millions)", fill = "Group"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(FIG_DIR, "fig_library_size.png"), p_lib,
       width = 8, height = 10, dpi = 150)
message("Saved: fig_library_size.png")

# ── 6. PCA ────────────────────────────────────────────────────────────────────
# Use top 1,000 most variable genes
gene_vars   <- apply(log2cpm_mat, 1, var, na.rm = TRUE)
top1k_genes <- names(sort(gene_vars, decreasing = TRUE))[1:1000]
pca_mat     <- t(log2cpm_mat[top1k_genes, ])

pca_res <- prcomp(pca_mat, center = TRUE, scale. = FALSE)
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$sample <- rownames(pca_df)
pca_df$group  <- col_data[pca_df$sample, "group"]

# Flag outliers (|PC1| > mean + 2*SD)
pc1_thresh <- mean(pca_df$PC1) + 2 * sd(pca_df$PC1)
pca_df$outlier <- abs(pca_df$PC1) > abs(pc1_thresh)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  stat_ellipse(aes(group = group), level = 0.80, linetype = "dashed", size = 0.6) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(
    data = subset(pca_df, outlier),
    aes(label = sample), size = 2.5, show.legend = FALSE
  ) +
  scale_color_manual(values = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  labs(
    title    = "PCA — GSE224615 PBMC RNA-seq (top 1,000 variable genes)",
    x        = paste0("PC1 (", pct_var[1], "% variance)"),
    y        = paste0("PC2 (", pct_var[2], "% variance)"),
    color    = "Group",
    caption  = "Dashed ellipses: 80% concentration regions"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig1_pca.png"), p_pca,
       width = 8, height = 6, dpi = 150)
message("Saved: fig1_pca.png")
message("  PC1: ", pct_var[1], "% | PC2: ", pct_var[2], "%")

# ── 7. Save processed objects ─────────────────────────────────────────────────
saveRDS(log2cpm_mat, file.path(RESULTS_DIR, "expr_mat.rds"))
saveRDS(col_data,    file.path(RESULTS_DIR, "col_data.rds"))
message("\nScript 01 complete.")
message("  expr_mat.rds  : ", nrow(log2cpm_mat), " genes x ",
        ncol(log2cpm_mat), " samples")
message("  col_data.rds  : ", nrow(col_data), " samples")
