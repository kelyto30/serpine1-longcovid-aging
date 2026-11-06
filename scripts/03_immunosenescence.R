# =============================================================================
# Script 03: Immunosenescence / exhaustion GSVA, marker genes, co-enrichment
# Project: Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with
#          NF-kB Activation in COVID-19 Postmortem Brain
# Author:  Kelechi Wisdom Elechi | PI: Gabriel A. de Erausquin
#          UT Health San Antonio / Barshop Institute
# -----------------------------------------------------------------------------
# PROVENANCE: The code below is the code that was executed to produce the
# results reported in the manuscript, taken from the R console log of that run
# (Phase 4 run). R 4.3.2, GSVA 1.50.
# Changes from the logged code (marked "# [v1.1.0]" where they occur):
#   - directory creation added
#   - NOTE: co-enrichment correlations use Pearson's r (method = 'pearson') across all 35 samples
# INPUT:  data/expr_outlier_removed.rds, data/gsva_senescence_scores*.rds
# OUTPUT: data/gsva_immunosenescence_scores*.rds, data/coenrichment_scores.rds; results/immunosenescence/*
# Run from the repository root.
# =============================================================================

for (d in c("data", "results/qc", "results/deg", "results/senescence",
            "results/immunosenescence", "results/integration")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)   # [v1.1.0] create output dirs
}

suppressPackageStartupMessages({
  library(GSVA); library(ggplot2); library(ggrepel)
  library(pheatmap); library(dplyr); library(tidyr)
  library(patchwork); library(RColorBrewer); library(ggbeeswarm)
})

cat("=== Loading objects ===\n")
expr_clean   <- readRDS("data/expr_outlier_removed.rds")
col_clean    <- readRDS("data/coldata_outlier_removed.rds")
sen_scores   <- readRDS("data/gsva_senescence_scores.rds")
scores_df    <- readRDS("data/gsva_senescence_scores_df.rds")

cat("Matrix:", nrow(expr_clean), "genes x", ncol(expr_clean), "samples\n")
cat("Groups:\n"); print(table(col_clean$group))

grp_cols  <- c("Non_LongCOVID" = "#2196F3", "LongCOVID" = "#F44336")
expressed <- rownames(expr_clean)

# ── STEP 1: DEFINE IMMUNOSENESCENCE GENE SETS ─────────────────────────────────
cat("\n=== Defining immunosenescence gene sets ===\n")

immunosen_panel <- c("KLRG1","B3GAT1","TIGIT","LAG3","PDCD1",
                     "HAVCR2","ENTPD1","FCGR3A","NCR1",
                     "SELL","CCR7","IL7R","TCF7")

exhaustion_tfs <- c("TOX","TOX2","NR4A1","NR4A2","NR4A3",
                    "BATF","IRF4","PRDM1","EOMES","TBX21")

exhaustion_sig <- c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4",
                    "TOX","TOX2","NR4A1","NR4A2","NR4A3",
                    "BATF","IRF4","PRDM1","EOMES","TBX21",
                    "ENTPD1","TCF7","GZMB","GZMK","PRF1",
                    "SELL","CCR7","IL7R","CD28","CD27",
                    "KLRG1","B3GAT1","CX3CR1","FCGR3A")

nk_immunosen <- c("FCGR3A","NCR1","KLRG1","B3GAT1","CX3CR1",
                   "FCER1G","TYROBP","LILRB1")

isen_sets <- list(
  Immunosenescence_Panel = intersect(unique(immunosen_panel), expressed),
  Exhaustion_TFs         = intersect(unique(exhaustion_tfs),  expressed),
  Exhaustion_Signature   = intersect(unique(exhaustion_sig),  expressed),
  NK_Immunosenescence    = intersect(unique(nk_immunosen),    expressed)
)

for (gs in names(isen_sets)) {
  cat(gs, ":", length(isen_sets[[gs]]), "genes in dataset\n")
}

# ── STEP 2: GSVA IMMUNOSENESCENCE SCORING ────────────────────────────────────
cat("\n=== Running GSVA immunosenescence scoring ===\n")

isen_param <- GSVA::gsvaParam(
  exprData = expr_clean,
  geneSets = isen_sets,
  minSize  = 5,
  maxSize  = 500
)

isen_scores <- GSVA::gsva(isen_param, verbose = FALSE)

cat("Immunosenescence scoring complete.\n")
print(round(isen_scores, 3))

isen_df <- as.data.frame(t(isen_scores))
isen_df$sample_id <- rownames(isen_df)
isen_df <- dplyr::left_join(isen_df, col_clean, by = "sample_id")

# ── STEP 3: DIRECT EXPRESSION — KEY IMMUNOSENESCENCE MARKERS ─────────────────
cat("\n=== Direct expression of immunosenescence markers ===\n")

key_isen_genes <- c("KLRG1","B3GAT1","TIGIT","LAG3","PDCD1",
                    "TOX","NR4A1","NR4A2","CD28","CD27",
                    "SELL","IL7R","TCF7","FCGR3A","GZMB")
key_isen_in <- intersect(key_isen_genes, expressed)
cat("Key immunosenescence genes in dataset:", length(key_isen_in), "\n")
cat("Genes found:", paste(key_isen_in, collapse=", "), "\n")

# Direct expression comparison for each marker
marker_stats <- lapply(key_isen_in, function(g) {
  lc  <- expr_clean[g, col_clean$sample_id[col_clean$group=="LongCOVID"]]
  nlc <- expr_clean[g, col_clean$sample_id[col_clean$group=="Non_LongCOVID"]]
  wt  <- wilcox.test(as.numeric(lc), as.numeric(nlc), exact=FALSE)
  data.frame(
    gene         = g,
    median_LC    = round(median(as.numeric(lc)),  3),
    median_NonLC = round(median(as.numeric(nlc)), 3),
    delta        = round(median(as.numeric(lc)) - median(as.numeric(nlc)), 3),
    p_value      = round(wt$p.value, 4)
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(padj = round(p.adjust(p_value, "BH"), 4)) %>%
  dplyr::arrange(p_value)

cat("\nImmunoscenescence marker gene statistics:\n")
print(marker_stats)

# ── STEP 4: STATISTICAL TESTING ON GSVA SCORES ───────────────────────────────
cat("\n=== GSVA score statistics: LongCOVID vs Non_LongCOVID ===\n")

isen_stats <- lapply(names(isen_sets), function(gs) {
  lc  <- isen_df[[gs]][isen_df$group == "LongCOVID"]
  nlc <- isen_df[[gs]][isen_df$group == "Non_LongCOVID"]
  wt  <- wilcox.test(lc, nlc, exact = FALSE)
  data.frame(
    gene_set     = gs,
    n_LC         = length(lc),
    n_NonLC      = length(nlc),
    median_LC    = round(median(lc),  4),
    median_NonLC = round(median(nlc), 4),
    delta        = round(median(lc) - median(nlc), 4),
    p_value      = round(wt$p.value, 4)
  )
}) %>% dplyr::bind_rows()
isen_stats$padj <- round(p.adjust(isen_stats$p_value, "BH"), 4)

cat("Immunosenescence GSVA results:\n")
print(isen_stats[, c("gene_set","delta","p_value","padj")])

# ── STEP 5: CO-ENRICHMENT TEST ────────────────────────────────────────────────
cat("\n=== KEY TEST: Immunosenescence vs Senescence co-enrichment ===\n")

coenrich_df <- data.frame(
  sample_id  = colnames(expr_clean),
  ImmunoSen  = as.numeric(isen_scores["Immunosenescence_Panel", ]),
  Exhaustion = as.numeric(isen_scores["Exhaustion_Signature",   ]),
  SenMayo    = as.numeric(sen_scores["SenMayo",  ]),
  SASP_Core  = as.numeric(sen_scores["SASP_Core",])
) %>% dplyr::left_join(col_clean, by = "sample_id")

cor_isen_sen <- cor.test(coenrich_df$ImmunoSen, coenrich_df$SenMayo,
                          method = "pearson")
cor_exh_sasp <- cor.test(coenrich_df$Exhaustion, coenrich_df$SASP_Core,
                          method = "pearson")
cor_isen_sasp <- cor.test(coenrich_df$ImmunoSen, coenrich_df$SASP_Core,
                           method = "pearson")

cat("Immunosenescence Panel vs SenMayo  : r =",
    round(cor_isen_sen$estimate, 3), " p =", round(cor_isen_sen$p.value, 4), "\n")
cat("Exhaustion Signature vs SASP_Core  : r =",
    round(cor_exh_sasp$estimate, 3), " p =", round(cor_exh_sasp$p.value, 4), "\n")
cat("Immunosenescence Panel vs SASP_Core: r =",
    round(cor_isen_sasp$estimate, 3), " p =", round(cor_isen_sasp$p.value, 4), "\n")

# Co-enrichment scatter plot
p_coenrich <- ggplot(coenrich_df,
                      aes(x=ImmunoSen, y=SenMayo, colour=group)) +
  geom_point(size=3.5, alpha=0.85) +
  geom_smooth(method="lm", se=TRUE, colour="black",
              linetype="dashed", linewidth=0.8) +
  ggrepel::geom_text_repel(aes(label=sample_id), size=2.2,
                             max.overlaps=15, segment.size=0.3) +
  scale_colour_manual(values=grp_cols, name="Group") +
  annotate("text",
           x = min(coenrich_df$ImmunoSen) + 0.05,
           y = max(coenrich_df$SenMayo)   - 0.05,
           label = paste0("r = ", round(cor_isen_sen$estimate, 3),
                           "\np = ", round(cor_isen_sen$p.value, 3)),
           size=4, fontface="bold", hjust=0) +
  labs(title    = "Co-Enrichment: Immunosenescence vs Cellular Senescence",
       subtitle  = "GSE224615 | Pearson correlation across all 35 samples",
       x = "Immunosenescence Panel Score (GSVA)",
       y = "SenMayo Senescence Score (GSVA)") +
  theme_classic(base_size=11) +
  theme(plot.title = element_text(face="bold"))

ggsave("results/immunosenescence/fig11_coenrichment.pdf",
       p_coenrich, width=9, height=6, dpi=300)
ggsave("results/immunosenescence/fig11_coenrichment.png",
       p_coenrich, width=9, height=6, dpi=300)
cat("Figure 11 saved.\n")

# ── STEP 6: IMMUNOSENESCENCE HEATMAP ─────────────────────────────────────────
ann_col <- data.frame(Group=col_clean$group, row.names=col_clean$sample_id)

if (length(key_isen_in) >= 3) {
  pheatmap::pheatmap(
    expr_clean[key_isen_in, ],
    scale             = "row",
    annotation_col    = ann_col,
    annotation_colors = list(Group = grp_cols),
    color = colorRampPalette(c("#2196F3","white","#F44336"))(100),
    cluster_rows = TRUE, cluster_cols = TRUE,
    show_rownames = TRUE, show_colnames = FALSE,
    fontsize_row  = 9,
    main     = "Immunosenescence Marker Gene Expression - GSE224615",
    filename = "results/immunosenescence/fig9_immunosenescence_heatmap.pdf",
    width=10, height=7
  )
  cat("Figure 9 (heatmap) saved.\n")
}

# ── STEP 7: SAVE ──────────────────────────────────────────────────────────────
saveRDS(isen_scores,   "data/gsva_immunosenescence_scores.rds")
saveRDS(isen_df,       "data/gsva_immunosenescence_scores_df.rds")
saveRDS(coenrich_df,   "data/coenrichment_scores.rds")

write.csv(isen_stats,   "results/immunosenescence/immunosenescence_statistics.csv",
          row.names=FALSE)
write.csv(marker_stats, "results/immunosenescence/marker_gene_statistics.csv",
          row.names=FALSE)
write.csv(coenrich_df,  "results/integration/coenrichment_scores.csv",
          row.names=FALSE)

cat("\n=================================================================\n")
cat("  PHASE 4 (Sub-Aim 2.2) COMPLETE\n")
cat("=================================================================\n")
cat("Immunosenescence GSVA results:\n")
print(isen_stats[, c("gene_set","delta","p_value","padj")])
cat("\nCo-enrichment correlations:\n")
cat("  ImmunoSen vs SenMayo  : r =", round(cor_isen_sen$estimate,3),
    "p =", round(cor_isen_sen$p.value,4), "\n")
cat("  ImmunoSen vs SASP_Core: r =", round(cor_isen_sasp$estimate,3),
    "p =", round(cor_isen_sasp$p.value,4), "\n")
cat("\nTop immunosenescence marker genes:\n")
print(head(marker_stats, 5))
cat("=================================================================\n")
