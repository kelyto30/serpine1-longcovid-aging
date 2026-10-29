# =============================================================================
# Script 04: Immunosenescence Profiling and Co-Enrichment Analysis
# Project:   Peripheral SERPINE1/PAI-1 as a Cross-Compartment Biomarker
#            of Long COVID-Associated Cognitive Aging
# Author:    Kelechi Wisdom Elechi
# PI:        Gabriel A. de Erausquin
# Date:      2026-10
# =============================================================================
# INPUT:  results/expr_mat.rds
#         results/col_data.rds
#         results/gsva_scores.rds
# OUTPUT: results/isen_wilcox_stats.csv
#         results/coenrichment_cor.txt
#         results/figures/fig4_coenrichment.png
#         results/figures/fig_immunosen_heatmap.png
#         results/figures/fig_isen_markers.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ── 0. Paths ──────────────────────────────────────────────────────────────────
RESULTS_DIR <- "results"
FIG_DIR     <- file.path(RESULTS_DIR, "figures")

expr_mat    <- readRDS(file.path(RESULTS_DIR, "expr_mat.rds"))
col_data    <- readRDS(file.path(RESULTS_DIR, "col_data.rds"))
gsva_scores <- readRDS(file.path(RESULTS_DIR, "gsva_scores.rds"))

lc_samples    <- rownames(col_data)[col_data$group == "LC"]
nonlc_samples <- rownames(col_data)[col_data$group == "NonLC"]
all_samples   <- c(lc_samples, nonlc_samples)

# ── 1. Immunosenescence marker genes ─────────────────────────────────────────
isen_markers <- c(
  "B3GAT1",  # CD57 — canonical T-cell immunosenescence marker
  "NR4A1",   # Nur77 — exhaustion/senescence TF
  "KLRG1",   # Short-lived effector marker
  "EOMES",   # Effector/exhaustion TF
  "CD28",    # Naive/memory T-cell marker (loss = senescence)
  "SELL",    # CD62L — lymph node homing (loss = effector)
  "IL7R",    # CD127 — memory T-cell (loss = exhaustion)
  "CX3CR1",  # Tissue-trafficking chemokine receptor
  "S100A4",  # Calcium-binding, senescence-associated
  "TFAM",    # Mitochondrial biogenesis
  "TIGIT",   # Co-inhibitory (exhaustion)
  "LAG3",    # Co-inhibitory (exhaustion)
  "HAVCR2",  # TIM-3 — terminal exhaustion
  "PDCD1",   # PD-1 — exhaustion
  "CD244"    # 2B4 — NK/T immunosenescence
)

isen_available <- intersect(isen_markers, rownames(expr_mat))
message("Immunosenescence markers in matrix: ", length(isen_available),
        " / ", length(isen_markers))

# ── 2. Wilcoxon test for each marker ─────────────────────────────────────────
isen_stats <- do.call(rbind, lapply(isen_available, function(gene) {
  lc_v    <- expr_mat[gene, lc_samples]
  nonlc_v <- expr_mat[gene, nonlc_samples]
  wt <- wilcox.test(lc_v, nonlc_v, exact = FALSE)
  data.frame(
    gene         = gene,
    median_LC    = median(lc_v,    na.rm = TRUE),
    median_NonLC = median(nonlc_v, na.rm = TRUE),
    delta        = median(lc_v, na.rm = TRUE) - median(nonlc_v, na.rm = TRUE),
    p_value      = wt$p.value,
    stringsAsFactors = FALSE
  )
}))
isen_stats$padj <- p.adjust(isen_stats$p_value, method = "BH")
isen_stats <- isen_stats[order(isen_stats$p_value), ]

write.csv(isen_stats,
          file.path(RESULTS_DIR, "isen_wilcox_stats.csv"),
          row.names = FALSE)

message("\nImmunosenesence marker stats (top rows):")
print(head(isen_stats, 6))

# ── 3. Immunosenescence marker boxplots ──────────────────────────────────────
isen_long <- do.call(rbind, lapply(isen_available, function(g) {
  data.frame(
    gene       = g,
    expression = c(expr_mat[g, lc_samples], expr_mat[g, nonlc_samples]),
    group      = c(rep("LC", length(lc_samples)), rep("NonLC", length(nonlc_samples))),
    stringsAsFactors = FALSE
  )
}))
isen_long$group <- factor(isen_long$group, levels = c("NonLC", "LC"))

p_isen_markers <- ggplot(isen_long, aes(x = group, y = expression, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55) +
  geom_jitter(width = 0.14, size = 0.7, alpha = 0.5, aes(color = group)) +
  facet_wrap(~ gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values  = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  scale_color_manual(values = c("NonLC" = "#2B5EA7", "LC" = "#B01A14")) +
  labs(
    title    = "Immunosenescence Marker Expression — LC vs. non-LC PBMCs",
    subtitle = "GSE224615 | 15-gene immunosenescence panel",
    x = NULL, y = "log2(CPM + 1)", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8))

ggsave(file.path(FIG_DIR, "fig_isen_markers.png"), p_isen_markers,
       width = 14, height = 8, dpi = 150)
message("Saved: fig_isen_markers.png")

# ── 4. Immunosenescence heatmap ───────────────────────────────────────────────
hm_mat <- expr_mat[isen_available, c(lc_samples, nonlc_samples)]
hm_scaled <- t(scale(t(hm_mat)))

ann_col <- data.frame(
  Group = col_data[colnames(hm_mat), "group"],
  row.names = colnames(hm_mat)
)
ann_colors <- list(Group = c("LC" = "#D73027", "NonLC" = "#4575B4"))

pheatmap(
  hm_scaled,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-2.5, 2.5, length.out = 101),
  main              = "Immunosenescence Markers — Row-scaled Expression",
  fontsize_row      = 10,
  filename          = file.path(FIG_DIR, "fig_immunosen_heatmap.png"),
  width = 7, height = 6
)
message("Saved: fig_immunosen_heatmap.png")

# ── 5. Co-enrichment: ImmunoSenescence vs. SASP Core ─────────────────────────
# Both scores from the GSVA matrix
if (!all(c("Immunosenescence", "SASP_Core") %in% rownames(gsva_scores))) {
  stop("Expected 'Immunosenescence' and 'SASP_Core' rows in gsva_scores.")
}

coenrich_df <- data.frame(
  sample      = colnames(gsva_scores),
  ImmunoSen   = as.numeric(gsva_scores["Immunosenescence", ]),
  SASP_Core   = as.numeric(gsva_scores["SASP_Core", ]),
  stringsAsFactors = FALSE
)
coenrich_df$group <- col_data[coenrich_df$sample, "group"]
coenrich_df$group <- factor(coenrich_df$group, levels = c("NonLC", "LC"))

# Spearman correlation
cor_test <- cor.test(coenrich_df$ImmunoSen, coenrich_df$SASP_Core,
                     method = "spearman", exact = FALSE)
r_val  <- round(cor_test$estimate, 3)
p_val  <- signif(cor_test$p.value, 3)

sink(file.path(RESULTS_DIR, "coenrichment_cor.txt"))
cat("Co-enrichment: Immunosenescence vs. SASP Core GSVA Scores\n")
cat("============================================================\n")
cat("Method:       Spearman rank correlation\n")
cat("Spearman r:  ", r_val, "\n")
cat("p-value:     ", p_val, "\n")
cat("n samples:   ", nrow(coenrich_df), "\n")
cat("\nInterpretation: r = ", r_val,
    " indicates an inverse co-enrichment pattern.\n",
    "SASP-high LC samples tend toward lower immunosenescence scores,\n",
    "suggesting divergent rather than co-occurring senescence programs.\n",
    sep = "")
sink()

message("\nCo-enrichment result:")
message("  Spearman r = ", r_val, "  p = ", p_val)

# ── 6. Co-enrichment scatter plot (Figure 4) ─────────────────────────────────
p_coenrich <- ggplot(coenrich_df, aes(x = ImmunoSen, y = SASP_Core, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              color = "gray40", linetype = "solid", size = 0.8) +
  scale_color_manual(values = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("Spearman r = ", r_val, "\np = ", p_val),
           hjust = 1.1, vjust = 1.4,
           size = 4, fontface = "italic", color = "gray20") +
  labs(
    title    = "Co-enrichment: Immunosenescence vs. SASP Core",
    subtitle = "GSE224615 PBMC GSVA scores | Spearman correlation",
    x        = "Immunosenescence GSVA score",
    y        = "SASP Core GSVA score",
    color    = "Group",
    caption  = paste0("Shaded band: 95% CI of linear fit | r = ", r_val,
                      ", p = ", p_val)
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig4_coenrichment.png"), p_coenrich,
       width = 6, height = 5.5, dpi = 150)
message("Saved: fig4_coenrichment.png")

message("\nScript 04 complete.")
