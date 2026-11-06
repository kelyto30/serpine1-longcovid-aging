# =============================================================================
# Script 02: Outlier removal + GSVA senescence scoring (SenMayo, CellAge, SASP Core, p21/p16/p53)
# Project: Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with
#          NF-kB Activation in COVID-19 Postmortem Brain
# Author:  Kelechi Wisdom Elechi
#          UT Health San Antonio / Barshop Institute
# -----------------------------------------------------------------------------
# PROVENANCE: The code below is the code that was executed to produce the
# results reported in the manuscript, taken from the R console log of that run
# (Phase 3 run; GSVA step re-run after installing GSVA, identical code). R 4.3.2, GSVA 1.50.
# Changes from the logged code (marked "# [v1.1.0]" where they occur):
#   - directory creation added
#   - NOTE: GSVA::gsvaParam() with default kcdf -> classic GSVA (Hanzelmann 2013), NOT ssGSEA
# INPUT:  data/expr_log2_symbol_gse224615.rds, data/coldata_gse224615.rds
# OUTPUT: data/gsva_senescence_scores*.rds, data/expr_outlier_removed.rds; results/senescence/*
# Run from the repository root.
# =============================================================================

for (d in c("data", "results/qc", "results/deg", "results/senescence",
            "results/immunosenescence", "results/integration")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)   # [v1.1.0] create output dirs
}


suppressPackageStartupMessages({
  library(GSVA); library(ggplot2); library(ggrepel)
  library(pheatmap); library(dplyr); library(tidyr)
  library(patchwork); library(RColorBrewer); library(viridis)
  library(ggbeeswarm); library(scales)
})

cat("=== Loading Phase 1 objects ===\n")
expr_sym     <- readRDS("data/expr_log2_symbol_gse224615.rds")
col_data     <- readRDS("data/coldata_gse224615.rds")
deg_results  <- readRDS("data/deseq2_results_LC_vs_NonLC.rds")

cat("Matrix:", nrow(expr_sym), "genes x", ncol(expr_sym), "samples\n")
cat("Groups:\n"); print(table(col_data$group))

grp_cols <- c("Non_LongCOVID" = "#2196F3", "LongCOVID" = "#F44336")

# ── STEP 1: REMOVE CONFIRMED OUTLIER ─────────────────────────────────────────
cat("\n=== Removing outlier GSM7027490 (PC1 dev = 4.42 SD) ===\n")

outlier_remove <- "GSM7027490"
expr_clean     <- expr_sym[, !colnames(expr_sym) %in% outlier_remove]
col_clean      <- col_data[!rownames(col_data) %in% outlier_remove, ]
col_clean$group <- droplevels(col_clean$group)

cat("Samples after outlier removal:", ncol(expr_clean), "\n")
cat("Groups:\n"); print(table(col_clean$group))

# ── STEP 2: DEFINE SENESCENCE GENE SETS ──────────────────────────────────────
cat("\n=== Defining senescence gene sets ===\n")

expressed <- rownames(expr_clean)

senmayo_genes <- c(
  "IL6","CXCL8","CXCL1","CXCL2","CXCL3","CXCL5","CCL2","CCL3","CCL4",
  "CCL5","CCL7","CCL8","CCL20","MMP1","MMP3","MMP9","MMP10","MMP13","MMP14",
  "IL1A","IL1B","IL7","IL13","IL15","IL33","CSF2","VEGFA","HGF","IGF1",
  "SERPINE1","SERPINE2","ICAM1","VCAM1","GDF15",
  "IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP6","IGFBP7",
  "CDKN1A","CDKN2A","CDKN2B","CDKN2D","TP53","RB1","CDKN1B","CDKN1C",
  "H2AX","BRCA1","ATM","ATR","CHEK1","CHEK2","MDM2","MDM4",
  "GLB1","LMNB1","HMGA1","HMGA2","BMI1","EZH2",
  "BCL2","BCL2L1","MCL1","SIRT1","SIRT3","SIRT6",
  "MTOR","AKT1","CGAS","STING1","IRF3",
  "WNT2","WNT5A","NOTCH1","TGM2","TGFB1","TGFB2","PTGS2"
)

cellage_genes <- c(
  "TP53","RB1","CDKN1A","CDKN2A","CDKN2B","CDKN2C","CDKN2D",
  "CDKN1B","CDKN1C","E2F1","E2F2","E2F3",
  "SIRT1","SIRT2","SIRT3","SIRT6","SIRT7",
  "LMNB1","LMNB2","LMNA","HMGB1","HMGB2","HMGA1","HMGA2",
  "GLB1","SERPINE1","IGFBP5","IGFBP7","GDF15",
  "IL6","CXCL8","MMP3","MMP9",
  "ATM","ATR","CHEK1","CHEK2",
  "BCL2","BCL2L1","BCL2L2","BAX","BBC3",
  "MDM2","MDM4","MTOR","AKT1","TERC","TERT","POT1"
)

sasp_core <- c("IL6","CXCL8","IL1A","IL1B","MMP3","MMP9",
                "CXCL1","CCL2","VEGFA","SERPINE1","ICAM1",
                "CSF2","GDF15","IGFBP7","IGFBP3","TNF")

p21_p16_p53 <- c("CDKN1A","CDKN2A","TP53","RB1","E2F1","MDM2","CDKN1B")

sen_sets <- list(
  SenMayo     = intersect(unique(senmayo_genes), expressed),
  CellAge     = intersect(unique(cellage_genes), expressed),
  SASP_Core   = intersect(unique(sasp_core),    expressed),
  p21_p16_p53 = intersect(unique(p21_p16_p53),  expressed)
)

for (gs in names(sen_sets)) {
  cat(gs, ":", length(sen_sets[[gs]]), "genes in dataset\n")
}

# ── STEP 3: GSVA SCORING ──────────────────────────────────────────────────────
cat("\n=== Running GSVA senescence scoring ===\n")

gsva_param <- GSVA::gsvaParam(
  exprData = expr_clean,
  geneSets = sen_sets,
  minSize  = 5,
  maxSize  = 500
)

gsva_scores <- GSVA::gsva(gsva_param, verbose = FALSE)

cat("GSVA complete. Score matrix:", nrow(gsva_scores), "sets x",
    ncol(gsva_scores), "samples\n")
print(round(gsva_scores, 3))

# Build scores dataframe
scores_df <- as.data.frame(t(gsva_scores))
scores_df$sample_id <- rownames(scores_df)
scores_df <- dplyr::left_join(scores_df, col_clean, by = "sample_id")

# ── STEP 4: STATISTICAL TESTING ───────────────────────────────────────────────
cat("\n=== Statistical comparison: LongCOVID vs Non_LongCOVID ===\n")

stats_results <- lapply(names(sen_sets), function(gs) {
  lc  <- scores_df[[gs]][scores_df$group == "LongCOVID"]
  nlc <- scores_df[[gs]][scores_df$group == "Non_LongCOVID"]
  wt  <- wilcox.test(lc, nlc, exact = FALSE)
  data.frame(
    gene_set       = gs,
    n_LC           = length(lc),
    n_NonLC        = length(nlc),
    median_LC      = round(median(lc),  4),
    median_NonLC   = round(median(nlc), 4),
    delta          = round(median(lc) - median(nlc), 4),
    W              = wt$statistic,
    p_value        = round(wt$p.value, 4)
  )
}) %>% dplyr::bind_rows()

stats_results$padj <- round(p.adjust(stats_results$p_value, method = "BH"), 4)

cat("\nSenescence Score Results:\n")
print(stats_results)

# ── STEP 5: BOXPLOTS ──────────────────────────────────────────────────────────
cat("\n=== Generating senescence score figures ===\n")

plot_score <- function(df, score_col, title, subtitle = NULL) {
  pval <- stats_results$p_value[stats_results$gene_set == score_col]
  padj <- stats_results$padj[stats_results$gene_set == score_col]
  plab <- paste0("p = ", round(pval, 3), " | padj = ", round(padj, 3))

  ggplot(df, aes(x = group, y = .data[[score_col]], fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.45, linewidth=0.5) +
    ggbeeswarm::geom_beeswarm(size = 2.8, alpha = 0.85, cex = 1.5) +
    scale_fill_manual(values = grp_cols, guide = "none") +
    annotate("text", x = 1.5, y = max(df[[score_col]]) * 1.05,
             label = plab, size = 3, fontface = "italic") +
    labs(title    = title,
         subtitle = subtitle,
         x = NULL, y = "GSVA Enrichment Score") +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 25, hjust = 1, size = 10),
          plot.title    = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8, colour = "gray40"))
}

p1 <- plot_score(scores_df, "SenMayo",
  "SenMayo Senescence Score",
  "105-gene validated senescence signature (Saul et al. 2022)")
p2 <- plot_score(scores_df, "CellAge",
  "CellAge Senescence Score",
  "Database-curated senescence gene set")
p3 <- plot_score(scores_df, "SASP_Core",
  "SASP Core Score",
  "Senescence-Associated Secretory Phenotype effectors")
p4 <- plot_score(scores_df, "p21_p16_p53",
  "p21/p16/p53 Cell Cycle Arrest Axis",
  "CDKN1A, CDKN2A, TP53, RB1, MDM2")

combined <- (p1 + p2) / (p3 + p4) +
  patchwork::plot_annotation(
    title    = "Sub-Aim 2.1: Senescence Transcriptomic Scores — GSE224615",
    subtitle = "Long COVID (n=22) vs Non-Long COVID (n=13) PBMCs | GSVA scoring | Wilcoxon test",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "gray40", size = 9)
    )
  )

ggsave("results/senescence/fig6_senescence_scores.pdf",
       combined, width = 12, height = 9, dpi = 300)
ggsave("results/senescence/fig6_senescence_scores.png",
       combined, width = 12, height = 9, dpi = 300)
cat("Figure 6 saved.\n")

# ── STEP 6: SENESCENCE GENE HEATMAP ──────────────────────────────────────────
cat("\n=== Generating senescence gene heatmap ===\n")

key_genes <- c("CDKN1A","CDKN2A","TP53","IL6","CXCL8","MMP3","MMP9",
               "SERPINE1","GDF15","IGFBP3","IGFBP7","GLB1","LMNB1",
               "HMGB1","BCL2L1","ATM","TGFB1","VEGFA","CCL2","IL1B")
key_genes_in <- intersect(key_genes, rownames(expr_clean))

heatmap_mat <- expr_clean[key_genes_in, ]
ann_col <- data.frame(Group = col_clean$group, row.names = col_clean$sample_id)

pheatmap::pheatmap(
  heatmap_mat,
  scale             = "row",
  annotation_col    = ann_col,
  annotation_colors = list(Group = grp_cols),
  color             = colorRampPalette(c("#2196F3","white","#F44336"))(100),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 9,
  main              = "Senescence Axis Gene Expression - GSE224615",
  filename          = "results/senescence/fig7_senescence_heatmap.pdf",
  width = 10, height = 7
)
cat("Figure 7 saved.\n")

# ── STEP 7: SERPINE1 SPOTLIGHT ────────────────────────────────────────────────
# SERPINE1 (PAI-1) was nominally significant in pre-computed DEGs (padj=0.282)
# AND is a canonical SASP gene — spotlight it directly
cat("\n=== SERPINE1 spotlight (key SASP gene) ===\n")

if ("SERPINE1" %in% rownames(expr_clean)) {
  serp_df <- data.frame(
    sample    = colnames(expr_clean),
    SERPINE1  = as.numeric(expr_clean["SERPINE1", ]),
    group     = col_clean$group
  )
  serp_wt <- wilcox.test(
    serp_df$SERPINE1[serp_df$group == "LongCOVID"],
    serp_df$SERPINE1[serp_df$group == "Non_LongCOVID"],
    exact = FALSE
  )
  cat("SERPINE1 (PAI-1) — LongCOVID vs Non_LongCOVID:\n")
  cat("  LongCOVID median:", round(median(serp_df$SERPINE1[serp_df$group=="LongCOVID"]), 3), "\n")
  cat("  Non_LongCOVID median:", round(median(serp_df$SERPINE1[serp_df$group=="Non_LongCOVID"]), 3), "\n")
  cat("  Wilcoxon p =", round(serp_wt$p.value, 4), "\n")
}

# ── STEP 8: SAVE ─────────────────────────────────────────────────────────────
cat("\n=== Saving Sub-Aim 2.1 results ===\n")

saveRDS(gsva_scores, "data/gsva_senescence_scores.rds")
saveRDS(scores_df,   "data/gsva_senescence_scores_df.rds")
saveRDS(expr_clean,  "data/expr_outlier_removed.rds")
saveRDS(col_clean,   "data/coldata_outlier_removed.rds")

write.csv(stats_results, "results/senescence/senescence_statistics.csv",
          row.names = FALSE)
write.csv(scores_df,     "results/senescence/gsva_scores_per_sample.csv",
          row.names = FALSE)

cat("\n=================================================================\n")
cat("  PHASE 3 (Sub-Aim 2.1) COMPLETE\n")
cat("=================================================================\n")
cat("Outlier removed  : GSM7027490\n")
cat("Final n          : LongCOVID=22, Non_LongCOVID=13\n")
cat("Gene sets scored : SenMayo, CellAge, SASP_Core, p21/p16/p53\n")
cat("\nSenescence score statistics:\n")
print(stats_results[, c("gene_set","delta","p_value","padj")])
cat("\nNext: Run immunosenescence scoring (Phase 4)\n")
cat("=================================================================\n")

