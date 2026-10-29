# =============================================================================
# Script 03: GSVA Senescence Scoring
# Project:   Peripheral SERPINE1/PAI-1 as a Cross-Compartment Biomarker
#            of Long COVID-Associated Cognitive Aging
# Author:    Kelechi Wisdom Elechi
# PI:        Gabriel A. de Erausquin
# Date:      2026-10
# =============================================================================
# INPUT:  results/expr_mat.rds
#         results/col_data.rds
# OUTPUT: results/gsva_scores.rds
#         results/gsva_stats.csv
#         results/figures/fig3_gsva_scores.png
#         results/figures/fig_senescence_heatmap.png
# NOTES:  Requires Bioconductor: BiocManager::install("GSVA")
# =============================================================================

suppressPackageStartupMessages({
  library(GSVA)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# ── 0. Paths ──────────────────────────────────────────────────────────────────
RESULTS_DIR <- "results"
FIG_DIR     <- file.path(RESULTS_DIR, "figures")

expr_mat <- readRDS(file.path(RESULTS_DIR, "expr_mat.rds"))
col_data <- readRDS(file.path(RESULTS_DIR, "col_data.rds"))

lc_samples    <- rownames(col_data)[col_data$group == "LC"]
nonlc_samples <- rownames(col_data)[col_data$group == "NonLC"]

# ── 1. Senescence Gene Sets ───────────────────────────────────────────────────

## SenMayo — Saul et al. 2022, Nature Communications
senmayo <- c(
  "ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2","BMP6",
  "C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26","CCL3",
  "CCL3L1","CCL4","CCL5","CCL7","CCL8","CCND1","CD55","CD9","CSF1",
  "CSF2","CSF2RB","CTSB","CXCL1","CXCL10","CXCL12","CXCL16","CXCL2",
  "CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR","EREG","ESM1",
  "ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMNN","HMGA1",
  "ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3","IGFBP4","IGFBP5",
  "IGFBP6","IGFBP7","IL10","IL13","IL15","IL18","IL1A","IL1B","IL2",
  "IL32","IL6","IL6ST","IL7","INHA","INHBA","ITGA2","ITPKA","JUN",
  "KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13","MMP14","MMP2",
  "MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF","PIGF","PLAT",
  "PLAU","PLAUR","PTGES","RPS6KA5","SCAMP4","SELPLG","SEMA3F","SERPINE1",
  "SERPINE2","SPP1","TFPI2","THBS1","THBS2","TNFRSF10C","TNFRSF11B",
  "TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA","VEGFC","VGF","WNT16","WNT2"
)

## CellAge — Avelar et al. 2020, Genome Biology
cellage <- c(
  "ABL1","ALDH1A3","ATF3","ATM","BARD1","BLM","BRCA1","BRCA2","BTG1",
  "CCNA2","CCND1","CCND2","CCNE1","CCNE2","CDK1","CDK2","CDK4","CDK6",
  "CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CHEK1","CHEK2",
  "CUL5","DAPK1","DAPK2","DCUN1D3","DEC1","E2F1","E2F3","EGR1","EPAS1",
  "FAS","FOXO1","FOXO3","GADD45A","GADD45B","HIC1","HIF1A","HIPK2",
  "HMGA1","ID1","ID2","ID3","JUN","LATS1","LATS2","MDM2","MDM4","MKI67",
  "MLH1","MSH2","MYC","MYCN","NBN","NOTCH1","NQO1","P14ARF","P21","PCNA",
  "PIGA","PLK1","PTEN","RAD51","RB1","RBBP8","RRAGA","SFN","SKP2","SIRT1",
  "SIRT2","SIRT3","SIRT6","SMAD2","SMAD3","SMAD4","SOD2","SP1","TP53",
  "TP63","TP73","TSC2","VHL","WRN","XPA","XRCC1","ZEB1"
)

## SASP Core — Coppé et al. 2010 + Basisty et al. 2020
sasp_core <- c(
  "SERPINE1","TNF","IL6","IL1B","CXCL8","MMP1","MMP3","MMP9","MMP10",
  "MMP12","CDKN1A","CDKN2A","TP53","VEGFA","HMGA1","CCL2","CCL5","CCL7",
  "CXCL1","CXCL2","CXCL10","ICAM1","VCAM1","IGFBP3","IGFBP5","IGFBP7",
  "FN1","COL1A1","GDF15","THBS1","PLAUR","TIMP1","TIMP2","IL10"
)

## p21/p16/p53 Pathway
p21_p16_p53 <- c(
  "CDKN1A","CDKN2A","TP53","MDM2","RB1","E2F1","GADD45A","GADD45B",
  "PCNA","CCND1","CDK4","CDK6","ATM","CHEK1","CHEK2","PLK1","MKI67",
  "FOXO3","SIRT1","SIRT6"
)

## Immunosenescence Panel
imm_sen <- c(
  "B3GAT1","NR4A1","KLRG1","EOMES","CD28","SELL","IL7R","CX3CR1",
  "S100A4","TFAM","TIGIT","LAG3","HAVCR2","PDCD1","CD244"
)

gene_sets <- list(
  SenMayo        = senmayo,
  CellAge        = cellage,
  SASP_Core      = sasp_core,
  p21_p16_p53    = p21_p16_p53,
  Immunosenescence = imm_sen
)

# Filter to genes present in expression matrix
gene_sets_filt <- lapply(gene_sets, function(gs) {
  keep <- intersect(gs, rownames(expr_mat))
  message("  ", deparse(substitute(gs)), ": ", length(keep), " / ", length(gs), " genes present")
  keep
})
# Re-report with names
for (nm in names(gene_sets_filt)) {
  message("  Gene set [", nm, "]: ", length(gene_sets_filt[[nm]]),
          " / ", length(gene_sets[[nm]]), " genes in matrix")
}

# ── 2. GSVA scoring ───────────────────────────────────────────────────────────
message("\nRunning GSVA (ssgsea method)...")
gsva_param <- gsvaParam(
  exprData   = as.matrix(expr_mat),
  geneSets   = gene_sets_filt,
  kcdf       = "Gaussian",   # appropriate for log2(CPM+1) continuous data
  minSize    = 5
)
gsva_scores <- gsva(gsva_param, verbose = FALSE)
message("GSVA complete. Score matrix: ", nrow(gsva_scores), " sets x ",
        ncol(gsva_scores), " samples")

saveRDS(gsva_scores, file.path(RESULTS_DIR, "gsva_scores.rds"))

# ── 3. Statistical testing ────────────────────────────────────────────────────
gsva_stats <- do.call(rbind, lapply(rownames(gsva_scores), function(gs) {
  lc_v    <- gsva_scores[gs, lc_samples]
  nonlc_v <- gsva_scores[gs, nonlc_samples]
  wt <- wilcox.test(lc_v, nonlc_v, exact = FALSE)
  data.frame(
    gene_set     = gs,
    mean_LC      = mean(lc_v),
    mean_NonLC   = mean(nonlc_v),
    delta        = mean(lc_v) - mean(nonlc_v),
    p_value      = wt$p.value,
    stringsAsFactors = FALSE
  )
}))
gsva_stats$padj <- p.adjust(gsva_stats$p_value, method = "BH")

write.csv(gsva_stats,
          file.path(RESULTS_DIR, "gsva_stats.csv"),
          row.names = FALSE)

message("\nGSVA Stats:")
print(gsva_stats[, c("gene_set", "delta", "p_value", "padj")])

# ── 4. Boxplot figure (Figure 3) ─────────────────────────────────────────────
gsva_long <- as.data.frame(t(gsva_scores)) %>%
  tibble::rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "gene_set", values_to = "score") %>%
  left_join(col_data %>% tibble::rownames_to_column("sample"),
            by = "sample")

# P-value label
p_labels <- setNames(
  paste0("p = ", signif(gsva_stats$p_value, 2)),
  gsva_stats$gene_set
)

gsva_long$group <- factor(gsva_long$group, levels = c("NonLC", "LC"))
gsva_long$p_label <- p_labels[gsva_long$gene_set]

p_gsva <- ggplot(gsva_long, aes(x = group, y = score, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.55) +
  geom_jitter(width = 0.14, size = 0.9, alpha = 0.55, aes(color = group)) +
  geom_text(
    data = gsva_long %>%
      group_by(gene_set, p_label) %>%
      summarise(y_pos = max(score) * 1.05, .groups = "drop") %>%
      distinct(),
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE, size = 3.2, hjust = 0.5
  ) +
  facet_wrap(~ gene_set, scales = "free_y", ncol = 3) +
  scale_fill_manual(values  = c("NonLC" = "#4575B4", "LC" = "#D73027")) +
  scale_color_manual(values = c("NonLC" = "#2B5EA7", "LC" = "#B01A14")) +
  labs(
    title    = "GSVA Senescence Enrichment Scores — LC vs. non-LC PBMCs",
    subtitle = "GSE224615 | Five curated senescence gene sets",
    x = NULL, y = "GSVA enrichment score", fill = NULL, color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text      = element_text(face = "bold", size = 9))

ggsave(file.path(FIG_DIR, "fig3_gsva_scores.png"), p_gsva,
       width = 12, height = 8, dpi = 150)
message("Saved: fig3_gsva_scores.png")

# ── 5. Heatmap of top senescence genes ───────────────────────────────────────
# Top 20 most variable senescence genes across all 5 sets
all_sen_genes <- unique(unlist(gene_sets_filt))
all_sen_genes <- intersect(all_sen_genes, rownames(expr_mat))
gene_var <- apply(expr_mat[all_sen_genes, ], 1, var, na.rm = TRUE)
top20 <- names(sort(gene_var, decreasing = TRUE))[1:20]

hm_mat <- expr_mat[top20, c(lc_samples, nonlc_samples)]
hm_mat_scaled <- t(scale(t(hm_mat)))  # row-scale (z-score)

ann_col <- data.frame(
  Group = col_data[colnames(hm_mat), "group"],
  row.names = colnames(hm_mat)
)
ann_colors <- list(Group = c("LC" = "#D73027", "NonLC" = "#4575B4"))

pheatmap(
  hm_mat_scaled,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks            = seq(-2.5, 2.5, length.out = 101),
  main              = "Top 20 Senescence Genes — Row-scaled Expression",
  fontsize_row      = 9,
  filename          = file.path(FIG_DIR, "fig_senescence_heatmap.png"),
  width = 8, height = 7
)
message("Saved: fig_senescence_heatmap.png")

message("\nScript 03 complete.")
