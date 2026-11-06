# =============================================================================
# Script 01: GEO metadata, expression matrix, QC, PCA, outlier detection
# Project: Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with
#          NF-kB Activation in COVID-19 Postmortem Brain
# Author:  Kelechi Wisdom Elechi
#          UT Health San Antonio / Barshop Institute
# -----------------------------------------------------------------------------
# PROVENANCE: The code below is the code that was executed to produce the
# results reported in the manuscript, taken from the R console log of that run
# (Phase 1 final run). R 4.3.2, GSVA 1.50.
# Changes from the logged code (marked "# [v1.1.0]" where they occur):
#   - directory creation added; 'Next' message points to renamed script
# INPUT:  data/GSE224615_DEGs.xlsx (GEO supplementary file); GEO series matrix via GEOquery
# OUTPUT: data/*.rds; results/qc/fig1-3; results/deg/precomputed_DEGs_LC_vs_NonLC.csv
# Run from the repository root.
# =============================================================================

for (d in c("data", "results/qc", "results/deg", "results/senescence",
            "results/immunosenescence", "results/integration")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)   # [v1.1.0] create output dirs
}


library(GEOquery); library(ggplot2); library(ggrepel)
library(pheatmap); library(dplyr); library(tidyr)
library(RColorBrewer); library(viridis); library(readxl)
library(org.Hs.eg.db); library(AnnotationDbi)

# ── RELOAD METADATA ───────────────────────────────────────────────────────────
gse     <- getGEO("GSE224615", GSEMatrix = TRUE, getGPL = FALSE)
gse_obj <- gse[[1]]
meta    <- Biobase::pData(gse_obj)

meta$group <- dplyr::case_when(
  trimws(meta[["lc status:ch1"]]) == "LC"     ~ "LongCOVID",
  trimws(meta[["lc status:ch1"]]) == "Non-LC" ~ "Non_LongCOVID",
  TRUE ~ NA_character_
)
meta$group <- factor(meta$group, levels = c("Non_LongCOVID","LongCOVID"))
meta$sex   <- dplyr::case_when(
  grepl("Female", meta[["characteristics_ch1.2"]]) ~ "F",
  grepl("Male",   meta[["characteristics_ch1.2"]]) ~ "M",
  TRUE ~ NA_character_
)
meta$pid <- trimws(gsub("pid: ", "", meta[["characteristics_ch1"]]))

# ── LOAD xlsx ─────────────────────────────────────────────────────────────────
xlsx_path <- "data/GSE224615_DEGs.xlsx"
raw_xl    <- readxl::read_excel(xlsx_path, sheet = 1)

# FIX 1: Exclude "Gene.name" from sample columns — it is a gene annotation column
stat_cols   <- c("ID","log2FoldChange","pvalue","padj","Gene.name")
sample_cols <- setdiff(colnames(raw_xl), stat_cols)
cat("Sample columns (n =", length(sample_cols), "):\n")
print(sample_cols)

# Pre-computed DEG table (use Gene.name for readable output)
deg_precomputed <- raw_xl %>%
  dplyr::select(ID, Gene.name, log2FoldChange, pvalue, padj) %>%
  dplyr::rename(ensembl = ID, gene = Gene.name) %>%
  dplyr::arrange(padj)

cat("Total genes:", nrow(deg_precomputed), "\n")
cat("Significant DEGs (padj<0.05):",
    sum(!is.na(deg_precomputed$padj) & deg_precomputed$padj < 0.05), "\n")
cat("Top 10 significant DEGs:\n")
print(head(deg_precomputed %>% dplyr::filter(!is.na(padj)), 10))

# ── BUILD EXPRESSION MATRIX ───────────────────────────────────────────────────
# FIX 2: Use Gene.name column as rownames (direct symbol mapping — no Ensembl conversion needed)
expr_mat_raw <- raw_xl %>%
  dplyr::select(all_of(c("Gene.name", sample_cols))) %>%
  dplyr::filter(!is.na(Gene.name), Gene.name != "") %>%
  dplyr::distinct(Gene.name, .keep_all = TRUE) %>%   # remove duplicate gene names
  tibble::column_to_rownames("Gene.name")

# FIX 3: Force numeric matrix
expr_mat <- apply(expr_mat_raw, 2, as.numeric)
rownames(expr_mat) <- rownames(expr_mat_raw)

cat("\nExpression matrix (gene symbols):", nrow(expr_mat), "genes x",
    ncol(expr_mat), "samples\n")
cat("Is numeric:", is.numeric(expr_mat), "\n")
cat("Any NAs:", sum(is.na(expr_mat)), "\n")

# ── MATCH SAMPLES TO METADATA ─────────────────────────────────────────────────
# FIX 4: Clean pid matching and remove Gene.name NA row
xlsx_pids_clean <- gsub("-.*", "", sample_cols)

pid_map <- data.frame(
  xlsx_col  = sample_cols,
  pid_clean = xlsx_pids_clean,
  stringsAsFactors = FALSE
)

meta_match <- meta %>%
  dplyr::select(geo_accession, pid, group, sex) %>%
  dplyr::inner_join(pid_map, by = c("pid" = "pid_clean"))  # inner_join drops unmatched

cat("\nMatched samples:", nrow(meta_match), "\n")
cat("Group distribution after matching:\n")
print(table(meta_match$group))

# Align expression matrix to matched metadata order
expr_aligned           <- expr_mat[, meta_match$xlsx_col]
colnames(expr_aligned) <- meta_match$geo_accession

# Build colData
col_data <- data.frame(
  sample_id = meta_match$geo_accession,
  pid       = meta_match$pid,
  group     = meta_match$group,
  sex       = meta_match$sex,
  row.names = meta_match$geo_accession,
  stringsAsFactors = FALSE
)
col_data$group <- factor(col_data$group,
                          levels = c("Non_LongCOVID","LongCOVID"))

cat("\ncolData:\n"); print(col_data)

# ── LOG2 TRANSFORM AND FILTER ─────────────────────────────────────────────────
cat("\n=== Preparing expression matrix for GSVA ===\n")

expr_log2 <- log2(expr_aligned + 1)

# Remove zero-variance genes
gene_var  <- apply(expr_log2, 1, var, na.rm = TRUE)
expr_filt <- expr_log2[gene_var > 0 & !is.na(gene_var), ]
cat("Genes after zero-variance filter:", nrow(expr_filt), "\n")

# Remove any remaining NA rows
expr_sym <- expr_filt[complete.cases(expr_filt), ]
cat("Genes after NA removal:", nrow(expr_sym), "\n")
cat("Final matrix:", nrow(expr_sym), "genes x", ncol(expr_sym), "samples\n")

# ── QC FIGURES ────────────────────────────────────────────────────────────────
cat("\n=== Generating QC figures ===\n")

grp_cols <- c("Non_LongCOVID" = "#2196F3", "LongCOVID" = "#F44336")

# Figure 1: Library size
lib_df <- data.frame(
  sample   = colnames(expr_aligned),
  lib_size = colSums(expr_aligned, na.rm = TRUE) / 1e6,
  group    = col_data$group
) %>% dplyr::arrange(group, lib_size)
lib_df$sample <- factor(lib_df$sample, levels = lib_df$sample)

p_lib <- ggplot(lib_df, aes(x = sample, y = lib_size, fill = group)) +
  geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = grp_cols, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title    = "Normalized Count Depth — GSE224615 (Long COVID PBMC)",
       subtitle = "Sum of normalized counts per sample (millions)",
       x = "Sample", y = "Count Sum (millions)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title  = element_text(face = "bold"))

ggsave("results/qc/fig1_library_sizes.pdf", p_lib, width=11, height=5, dpi=300)
ggsave("results/qc/fig1_library_sizes.png", p_lib, width=11, height=5, dpi=300)
cat("Figure 1 saved.\n")

# Figure 2: PCA
pca_res <- prcomp(t(expr_sym), scale. = FALSE)
pca_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_df <- data.frame(
  PC1    = pca_res$x[,1],
  PC2    = pca_res$x[,2],
  group  = col_data$group,
  sample = colnames(expr_sym)
)

p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, colour=group)) +
  geom_point(size=4, alpha=0.85) +
  ggrepel::geom_text_repel(aes(label=sample), size=2.5,
                             max.overlaps=20, segment.size=0.3) +
  stat_ellipse(level=0.80, linetype="dashed", linewidth=0.5) +
  scale_colour_manual(values=grp_cols, name="Group") +
  labs(title   = "PCA — GSE224615 Long COVID PBMC RNA-seq",
       subtitle = "log2(normalized counts + 1) | 80% ellipses",
       x = paste0("PC1 (", pca_var[1], "% variance)"),
       y = paste0("PC2 (", pca_var[2], "% variance)")) +
  theme_classic(base_size=11) +
  theme(plot.title = element_text(face="bold"))

ggsave("results/qc/fig2_pca.pdf", p_pca, width=9, height=6, dpi=300)
ggsave("results/qc/fig2_pca.png", p_pca, width=9, height=6, dpi=300)
cat("Figure 2 saved.\n")

# Figure 3: Sample distance heatmap
samp_dist     <- dist(t(expr_sym))
samp_dist_mat <- as.matrix(samp_dist)
ann_df        <- data.frame(Group = col_data$group,
                             row.names = col_data$sample_id)

pheatmap::pheatmap(
  samp_dist_mat,
  clustering_distance_rows  = samp_dist,
  clustering_distance_cols  = samp_dist,
  annotation_col            = ann_df,
  annotation_row            = ann_df,
  annotation_colors         = list(Group = grp_cols),
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")))(100),
  show_rownames  = FALSE, show_colnames = FALSE,
  main           = "Sample-to-Sample Distance — GSE224615",
  filename       = "results/qc/fig3_sample_distance_heatmap.pdf",
  width=8, height=7
)
cat("Figure 3 saved.\n")

# ── OUTLIER DETECTION ─────────────────────────────────────────────────────────
cat("\n=== Outlier detection ===\n")

outlier_df <- pca_df %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    pc1_dev = abs(PC1 - mean(PC1)) / sd(PC1),
    pc2_dev = abs(PC2 - mean(PC2)) / sd(PC2),
    outlier = (pc1_dev > 3 | pc2_dev > 3)
  ) %>% dplyr::ungroup()

outliers <- outlier_df %>% dplyr::filter(outlier)
if (nrow(outliers) == 0) {
  cat("No outliers detected — all samples within 3 SD of group centroid.\n")
} else {
  cat("Outliers detected:\n"); print(outliers[, c("sample","group","pc1_dev","pc2_dev")])
}

# ── SAVE ALL OBJECTS ──────────────────────────────────────────────────────────
cat("\n=== Saving objects ===\n")

saveRDS(expr_sym,        "data/expr_log2_symbol_gse224615.rds")
saveRDS(expr_aligned,    "data/expr_normalized_gse224615.rds")
saveRDS(col_data,        "data/coldata_gse224615.rds")
saveRDS(meta_match,      "data/metadata_gse224615.rds")
saveRDS(deg_precomputed, "data/deseq2_results_LC_vs_NonLC.rds")

write.csv(col_data,        "results/qc/sample_metadata_cleaned.csv")
write.csv(deg_precomputed, "results/deg/precomputed_DEGs_LC_vs_NonLC.csv",
          row.names = FALSE)

cat("\n=================================================================\n")
cat("  PHASE 1 COMPLETE\n")
cat("=================================================================\n")
cat("Dataset          : GSE224615\n")
cat("Samples          :", ncol(expr_sym), "\n")
cat("  Non_LongCOVID  :", sum(col_data$group == "Non_LongCOVID"), "\n")
cat("  LongCOVID      :", sum(col_data$group == "LongCOVID"), "\n")
cat("Genes (symbols)  :", nrow(expr_sym), "\n")
cat("Pre-computed DEGs (padj<0.05):",
    sum(!is.na(deg_precomputed$padj) & deg_precomputed$padj < 0.05), "\n")
cat("Outliers flagged :", nrow(outliers), "\n")
cat("Figures → results/qc/\n")
cat("Objects → data/\n")
cat("\nNext: Run scripts/02_senescence_scoring.R\n")  # [v1.1.0] script renamed
cat("=================================================================\n")

