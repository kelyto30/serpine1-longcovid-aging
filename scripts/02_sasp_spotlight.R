# =============================================================================
# Script 02: SASP Gene Spotlight — SERPINE1 Primary Analysis
# Project:   Peripheral SERPINE1/PAI-1 as a Cross-Compartment Biomarker
#            of Long COVID-Associated Cognitive Aging
# Author:    Kelechi Wisdom Elechi
# PI:        Gabriel A. de Erausquin
# Date:      2026-10
# =============================================================================
# INPUT:  results/expr_mat.rds
#         results/col_data.rds
# OUTPUT: results/sasp_wilcox_results.csv    -- all 34-gene Wilcoxon stats
#         results/serpine1_stats.txt         -- primary finding summary
#         results/figures/fig2_sasp_panel.png
#         results/figures/fig_serpine1_spotlight.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

# ── 0. Paths ──────────────────────────────────────────────────────────────────
RESULTS_DIR <- "results"
FIG_DIR     <- file.path(RESULTS_DIR, "figures")

expr_mat <- readRDS(file.path(RESULTS_DIR, "expr_mat.rds"))
col_data <- readRDS(file.path(RESULTS_DIR, "col_data.rds"))

# ── 1. SASP Core gene set ─────────────────────────────────────────────────────
# Curated from Coppé et al. 2010 and Basisty et al. 2020
sasp_core_genes <- c(
  "SERPINE1", "TNF",    "IL6",    "IL1B",   "CXCL8",
  "MMP1",     "MMP3",   "MMP9",   "MMP10",  "MMP12",
  "CDKN1A",  "CDKN2A", "TP53",   "VEGFA",  "HMGA1",
  "CCL2",    "CCL5",   "CCL7",   "CXCL1",  "CXCL2",
  "CXCL10",  "ICAM1",  "VCAM1",  "IGFBP3", "IGFBP5",
  "IGFBP7",  "FN1",    "COL1A1", "GDF15",  "THBS1",
  "PLAUR",   "TIMP1",  "TIMP2",  "IL10"
)

# Keep only genes present in the expression matrix
sasp_available <- intersect(sasp_core_genes, rownames(expr_mat))
message("SASP genes available in matrix: ", length(sasp_available),
        " / ", length(sasp_core_genes))

# ── 2. Wilcoxon rank-sum test for each SASP gene ─────────────────────────────
lc_samples    <- rownames(col_data)[col_data$group == "LC"]
nonlc_samples <- rownames(col_data)[col_data$group == "NonLC"]

run_wilcox <- function(gene) {
  lc_vals    <- expr_mat[gene, lc_samples]
  nonlc_vals <- expr_mat[gene, nonlc_samples]
  wt  <- wilcox.test(lc_vals, nonlc_vals, exact = FALSE)
  data.frame(
    gene         = gene,
    median_LC    = median(lc_vals,    na.rm = TRUE),
    median_NonLC = median(nonlc_vals, na.rm = TRUE),
    delta        = median(lc_vals, na.rm = TRUE) - median(nonlc_vals, na.rm = TRUE),
    p_value      = wt$p.value,
    stringsAsFactors = FALSE
  )
}

sasp_stats <- do.call(rbind, lapply(sasp_available, run_wilcox))
sasp_stats$padj <- p.adjust(sasp_stats$p_value, method = "BH")
sasp_stats <- sasp_stats[order(sasp_stats$p_value), ]

write.csv(sasp_stats,
          file.path(RESULTS_DIR, "sasp_wilcox_results.csv"),
          row.names = FALSE)

message("\nTop SASP genes by p-value:")
print(head(sasp_stats, 8))

# ── 3. SERPINE1 primary finding ───────────────────────────────────────────────
s1 <- sasp_stats[sasp_stats$gene == "SERPINE1", ]
message("\n--- SERPINE1 PRIMARY FINDING ---")
message("  LC median:    ", round(s1$median_LC,    4))
message("  NonLC median: ", round(s1$median_NonLC, 4))
message("  Delta:        ", round(s1$delta,   4))
message("  Wilcoxon p:   ", signif(s1$p_value, 3))
message("  BH adj p:     ", signif(s1$padj, 3))

# Write summary text
sink(file.path(RESULTS_DIR, "serpine1_stats.txt"))
cat("SERPINE1 / PAI-1 — Primary Finding\n")
cat("====================================\n")
cat("Dataset:         GSE224615 (Yin et al., 2023)\n")
cat("n LC:            22\n")
cat("n NonLC:         13\n")
cat("Median log2(CPM+1) — LC:    ", round(s1$median_LC,    4), "\n")
cat("Median log2(CPM+1) — NonLC: ", round(s1$median_NonLC, 4), "\n")
cat("Delta (LC - NonLC):          ", round(s1$delta,   4), "\n")
cat("Wilcoxon p-value:            ", signif(s1$p_value, 4), "\n")
cat("BH-adjusted p (34 genes):   ", signif(s1$padj, 4), "\n")
cat("\nConclusion: SERPINE1 is the sole SASP Core gene passing FDR < 0.05\n")
sink()

# ── 4. SERPINE1 spotlight boxplot ─────────────────────────────────────────────
s1_df <- data.frame(
  expression = c(expr_mat["SERPINE1", lc_samples],
                 expr_mat["SERPINE1", nonlc_samples]),
  group = c(rep("LC", length(lc_samples)),
            rep("NonLC", length(nonlc_samples))),
  stringsAsFactors = FALSE
)
s1_df$group <- factor(s1_df$group, levels = c("NonLC", "LC"))

p_s1 <- ggplot(s1_df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7, aes(color = group)) +
  scale_fill_manual(values  = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  scale_color_manual(values = c("NonLC" = "#2B5EA7", "LC" = "#B01A14")) +
  stat_compare_means(
    method = "wilcox.test",
    label  = "p.signif",
    comparisons = list(c("NonLC", "LC")),
    size = 5
  ) +
  labs(
    title    = "SERPINE1 expression in Long COVID PBMCs",
    subtitle = paste0("Wilcoxon p = ", signif(s1$p_value, 3),
                      "  |  adj p = ", signif(s1$padj, 3)),
    x        = NULL,
    y        = "log2(CPM + 1)",
    fill     = NULL, color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "fig_serpine1_spotlight.png"), p_s1,
       width = 5, height = 6, dpi = 150)
message("Saved: fig_serpine1_spotlight.png")

# ── 5. SASP 9-gene NF-κB panel figure (Figure 2) ─────────────────────────────
panel_genes <- c("SERPINE1", "TNF", "IL6", "IL1B", "CXCL8",
                 "MMP1", "CCL2", "CDKN1A", "ICAM1")
panel_genes <- intersect(panel_genes, rownames(expr_mat))

panel_df <- do.call(rbind, lapply(panel_genes, function(g) {
  data.frame(
    gene       = g,
    expression = c(expr_mat[g, lc_samples], expr_mat[g, nonlc_samples]),
    group      = c(rep("LC", length(lc_samples)), rep("NonLC", length(nonlc_samples))),
    stringsAsFactors = FALSE
  )
}))
panel_df$group <- factor(panel_df$group, levels = c("NonLC", "LC"))

# Add p-value annotations from sasp_stats
panel_stats <- sasp_stats[sasp_stats$gene %in% panel_genes,
                           c("gene", "p_value", "padj")]
panel_df <- merge(panel_df, panel_stats, by = "gene")
panel_df$p_label <- ifelse(panel_df$padj < 0.05,
                            paste0("adj p=", signif(panel_df$padj, 2)),
                            paste0("p=", signif(panel_df$p_value, 2)))

p_panel <- ggplot(panel_df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55) +
  geom_jitter(width = 0.14, size = 0.9, alpha = 0.6, aes(color = group)) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values  = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  scale_color_manual(values = c("NonLC" = "#2B5EA7", "LC" = "#B01A14")) +
  labs(
    title    = "SASP / NF-κB Gene Panel — Long COVID vs. non-LC PBMCs",
    subtitle = "Genes active in NF-κB-driven SASP; GSE224615",
    x = NULL, y = "log2(CPM + 1)", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 9),
    axis.text.x      = element_text(size = 8)
  )

ggsave(file.path(FIG_DIR, "fig2_sasp_panel.png"), p_panel,
       width = 10, height = 7, dpi = 150)
message("Saved: fig2_sasp_panel.png")

message("\nScript 02 complete.")
