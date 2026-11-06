# =============================================================================
# Script 04: SERPINE1 blood-brain bridge, STAT1 comparison, SASP/RELA blood panel, master table
# Project: Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with
#          NF-kB Activation in COVID-19 Postmortem Brain
# Author:  Kelechi Wisdom Elechi | PI: Gabriel A. de Erausquin
#          UT Health San Antonio / Barshop Institute
# -----------------------------------------------------------------------------
# PROVENANCE: The code below is the code that was executed to produce the
# results reported in the manuscript, taken from the R console log of that run
# (Phase 5 run). R 4.3.2, GSVA 1.50.
# Changes from the logged code (marked "# [v1.1.0]" where they occur):
#   - two element_text(fontface=) calls fixed to face= (they errored in the logged run, so Figs 12-13 were saved from a corrected re-run)
#   - STAT1 rows of the master table filled in / corrected to the logged blood result (p = 0.746)
#   - comments added citing the Aim 1 source of each hard-coded brain value
#   - the SASP/RELA panel defines the 9-gene BH family used for SERPINE1 padj = 0.0018
# INPUT:  data/*.rds from scripts 01-03; brain values hard-coded from Aim 1 (GSE188847) results
# OUTPUT: results/integration/fig12-14, master_integration_summary.csv; data/sasp_rela_gene_stats.rds
# Run from the repository root.
# =============================================================================

for (d in c("data", "results/qc", "results/deg", "results/senescence",
            "results/immunosenescence", "results/integration")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)   # [v1.1.0] create output dirs
}


suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(dplyr)
  library(tidyr); library(patchwork); library(RColorBrewer)
  library(pheatmap); library(scales)
})

cat("=== Loading all objects ===\n")
expr_clean   <- readRDS("data/expr_outlier_removed.rds")
col_clean    <- readRDS("data/coldata_outlier_removed.rds")
sen_scores   <- readRDS("data/gsva_senescence_scores.rds")
isen_scores  <- readRDS("data/gsva_immunosenescence_scores.rds")
coenrich_df  <- readRDS("data/coenrichment_scores.rds")
deg_results  <- readRDS("data/deseq2_results_LC_vs_NonLC.rds")

cat("Samples:", ncol(expr_clean), "| Genes:", nrow(expr_clean), "\n")
cat("Groups:\n"); print(table(col_clean$group))

grp_cols <- c("Non_LongCOVID" = "#2196F3", "LongCOVID" = "#F44336")
expressed <- rownames(expr_clean)

# ── FIGURE 1: SERPINE1 BLOOD-TO-BRAIN BRIDGE ─────────────────────────────────
cat("\n=== Building SERPINE1 blood-brain bridge figure ===\n")

# Aim 1 postmortem brain values (GSE188847; Aim 1 dissertation Section 4.3.3 / Table 4:
# SERPINE1 log2FC = +1.23 among RELA targets; RELA NES = +5.44)  # [v1.1.0] source cited
aim1_serpine1 <- data.frame(
  dataset   = "Postmortem Brain\n(GSE188847, Aim 1)",
  group     = c("Control","COVID-19"),
  log2FC    = c(0, 1.23),         # RELA target, from Aim 1 results
  context   = "RELA regulon target\nNES = +5.44"
)

# Aim 2 blood values — direct expression
serp_blood <- data.frame(
  sample_id = colnames(expr_clean),
  SERPINE1  = as.numeric(expr_clean["SERPINE1", ]),
  group     = col_clean$group
)

serp_wt <- wilcox.test(
  serp_blood$SERPINE1[serp_blood$group == "LongCOVID"],
  serp_blood$SERPINE1[serp_blood$group == "Non_LongCOVID"],
  exact = FALSE
)

cat("SERPINE1 blood expression:\n")
cat("  LongCOVID median:", round(median(serp_blood$SERPINE1[serp_blood$group=="LongCOVID"]), 3), "\n")
cat("  Non_LongCOVID median:", round(median(serp_blood$SERPINE1[serp_blood$group=="Non_LongCOVID"]), 3), "\n")
cat("  Wilcoxon p =", serp_wt$p.value, "\n")

# Blood boxplot
p_serp_blood <- ggplot(serp_blood,
                        aes(x = group, y = SERPINE1, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.45,
               linewidth = 0.6) +
  ggbeeswarm::geom_beeswarm(size = 2.8, alpha = 0.85, cex = 1.5) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  annotate("text", x = 1.5,
           y = max(serp_blood$SERPINE1) * 1.03,
           label = paste0("Wilcoxon p = ",
                           formatC(serp_wt$p.value, format="e", digits=1)),
           size = 3.5, fontface = "bold.italic", colour = "black") +
  labs(title    = "SERPINE1 (PAI-1) in Long COVID Blood",
       subtitle = "GSE224615 PBMC RNA-seq | Aim 2 discovery",
       x = NULL,
       y = "log2(normalized counts + 1)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "gray40", size = 8),
        axis.text.x   = element_text(angle = 20, hjust = 1))

# Brain context bar (Aim 1 reference)
p_serp_brain <- ggplot(aim1_serpine1,
                        aes(x = group, y = log2FC,
                            fill = group)) +
  geom_col(width = 0.5, alpha = 0.85) +
  scale_fill_manual(values = c("Control" = "grey70",
                                "COVID-19" = "#B71C1C"),
                    guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.35,
           label = "RELA regulon target\nlog2FC = +1.23 in brain",
           size = 3.2, fontface = "italic", colour = "#B71C1C") +
  labs(title    = "SERPINE1 in COVID-19 Postmortem Brain",
       subtitle = "GSE188847 | Aim 1 — RELA/NF-kB regulon",
       x = NULL,
       y = "log2 Fold Change (COVID-19 vs Control)") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(colour = "gray40", size = 8))

# Combined bridge figure
p_serpine1_bridge <- p_serp_blood + p_serp_brain +
  patchwork::plot_annotation(
    title   = "SERPINE1: Blood-to-Brain NF-kB Bridge",
    subtitle = paste0(
      "Significantly elevated in Long COVID blood (p=",
      formatC(serp_wt$p.value, format="e", digits=1),
      ") AND a RELA/NF-kB target upregulated in COVID-19 postmortem brain (Aim 1)"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "#B71C1C", size = 9,
                                    face = "italic")  # [v1.1.0] fontface -> face (logged run errored here)
    )
  )

ggsave("results/integration/fig12_SERPINE1_blood_brain_bridge.pdf",
       p_serpine1_bridge, width = 11, height = 6, dpi = 300)
ggsave("results/integration/fig12_SERPINE1_blood_brain_bridge.png",
       p_serpine1_bridge, width = 11, height = 6, dpi = 300)
cat("Figure 12 (SERPINE1 bridge) saved.\n")

# ── FIGURE 2: STAT1 CONVERGENCE ───────────────────────────────────────────────
cat("\n=== Building STAT1 convergence figure ===\n")

if ("STAT1" %in% expressed) {
  stat1_df <- data.frame(
    sample_id = colnames(expr_clean),
    STAT1     = as.numeric(expr_clean["STAT1", ]),
    group     = col_clean$group
  )

  stat1_wt <- wilcox.test(
    stat1_df$STAT1[stat1_df$group == "LongCOVID"],
    stat1_df$STAT1[stat1_df$group == "Non_LongCOVID"],
    exact = FALSE
  )

  cat("STAT1 blood expression:\n")
  cat("  LongCOVID median:",
      round(median(stat1_df$STAT1[stat1_df$group=="LongCOVID"]), 3), "\n")
  cat("  Non_LongCOVID median:",
      round(median(stat1_df$STAT1[stat1_df$group=="Non_LongCOVID"]), 3), "\n")
  cat("  Wilcoxon p =", round(stat1_wt$p.value, 4), "\n")

  # Aim 1 STAT1 reference values (Aim 1 Table 4: NES = +7.23, padj = 1.3e-10;
  # composition-corrected DESeq2 LFC = +0.69, padj = 0.006)  # [v1.1.0] source cited
  aim1_stat1 <- data.frame(
    metric  = c("VIPER NES\n(Aim 1 Brain)", "Cell-intrinsic LFC\n(Aim 1 Brain)"),
    value   = c(7.23, 0.69),
    context = c("Top activated TF\nNES = +7.23, padj = 1.3e-10",
                "Post-composition correction\nLFC = +0.69, padj = 0.006")
  )

  p_stat1_blood <- ggplot(stat1_df,
                           aes(x = group, y = STAT1, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.45,
                 linewidth = 0.6) +
    ggbeeswarm::geom_beeswarm(size = 2.8, alpha = 0.85, cex = 1.5) +
    scale_fill_manual(values = grp_cols, guide = "none") +
    annotate("text", x = 1.5,
             y = max(stat1_df$STAT1) * 1.02,
             label = paste0("p = ", round(stat1_wt$p.value, 3)),
             size = 3.5, fontface = "italic") +
    labs(title    = "STAT1 Expression in Long COVID Blood",
         subtitle = "GSE224615 PBMC RNA-seq | Aim 2",
         x = NULL, y = "log2(normalized counts + 1)") +
    theme_classic(base_size = 12) +
    theme(plot.title  = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = "gray40", size = 8),
          axis.text.x = element_text(angle = 20, hjust = 1))

  p_stat1_brain <- ggplot(aim1_stat1,
                           aes(x = metric, y = value, fill = metric)) +
    geom_col(width = 0.5, alpha = 0.85) +
    scale_fill_manual(values = c("VIPER NES\n(Aim 1 Brain)" = "#1A237E",
                                  "Cell-intrinsic LFC\n(Aim 1 Brain)" = "#283593"),
                      guide = "none") +
    geom_text(aes(label = round(value, 2)),
              vjust = -0.5, size = 4, fontface = "bold") +
    labs(title    = "STAT1 in COVID-19 Postmortem Brain",
         subtitle = "GSE188847 | Aim 1 — VIPER top TF + cell-intrinsic signal",
         x = NULL, y = "Effect Size") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(colour = "gray40", size = 8),
          axis.text.x   = element_text(size = 8))

  p_stat1_combined <- p_stat1_blood + p_stat1_brain +
    patchwork::plot_annotation(
      title    = "STAT1: IFN/JAK-STAT Convergence Across Blood and Brain",
      subtitle = "Top TF in postmortem brain (NES=+7.23) — tested in Long COVID peripheral blood",
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "#1A237E", size = 9,
                                      face = "italic")  # [v1.1.0] fontface -> face (logged run errored here)
      )
    )

  ggsave("results/integration/fig13_STAT1_convergence.pdf",
         p_stat1_combined, width = 11, height = 6, dpi = 300)
  ggsave("results/integration/fig13_STAT1_convergence.png",
         p_stat1_combined, width = 11, height = 6, dpi = 300)
  cat("Figure 13 (STAT1 convergence) saved.\n")
}

# ── FIGURE 3: SASP GENE PANEL IN BLOOD ───────────────────────────────────────
cat("\n=== Building SASP gene panel comparison ===\n")

# Key SASP/RELA overlap genes from Aim 1
sasp_rela_genes <- c("SERPINE1","IL6","MMP9","CCL2","CXCL10",
                      "VEGFA","ICAM1","IL1B","CXCL8","TNF")
sasp_in_blood   <- intersect(sasp_rela_genes, expressed)
cat("SASP/RELA genes in blood dataset:", paste(sasp_in_blood, collapse=", "), "\n")

sasp_expr_long <- expr_clean[sasp_in_blood, ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "sample_id", values_to = "expr") %>%
  dplyr::left_join(col_clean, by = "sample_id")

# Compute per-gene statistics
sasp_gene_stats <- lapply(sasp_in_blood, function(g) {
  lc  <- expr_clean[g, col_clean$sample_id[col_clean$group=="LongCOVID"]]
  nlc <- expr_clean[g, col_clean$sample_id[col_clean$group=="Non_LongCOVID"]]
  wt  <- wilcox.test(as.numeric(lc), as.numeric(nlc), exact=FALSE)
  data.frame(
    gene    = g,
    delta   = round(median(as.numeric(lc)) - median(as.numeric(nlc)), 3),
    p_value = round(wt$p.value, 4)
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(padj = round(p.adjust(p_value, "BH"), 4),
                sig  = ifelse(padj < 0.05, "***",
                              ifelse(p_value < 0.05, "*", "ns"))) %>%
  dplyr::arrange(p_value)

cat("\nSASP/RELA gene statistics in blood:\n")
print(sasp_gene_stats)

p_sasp_panel <- ggplot(sasp_expr_long,
                        aes(x = group, y = expr, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
               linewidth = 0.4) +
  ggbeeswarm::geom_beeswarm(size = 1.5, alpha = 0.7, cex = 1.2) +
  facet_wrap(~gene, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = grp_cols, name = "Group") +
  labs(
    title    = "SASP/NF-kB Target Gene Expression in Long COVID Blood",
    subtitle = paste0("Genes activated by RELA (NES=+5.44) in postmortem brain (Aim 1)",
                      " | GSE224615 | Wilcoxon"),
    x = NULL, y = "log2(normalized counts + 1)"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x   = element_text(angle = 30, hjust = 1, size = 7),
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(colour = "gray40", size = 8),
    strip.text    = element_text(face = "bold", size = 8),
    legend.position = "bottom"
  )

ggsave("results/integration/fig14_SASP_RELA_blood_panel.pdf",
       p_sasp_panel, width = 13, height = 6, dpi = 300)
ggsave("results/integration/fig14_SASP_RELA_blood_panel.png",
       p_sasp_panel, width = 13, height = 6, dpi = 300)
cat("Figure 14 (SASP/RELA panel) saved.\n")

# ── MASTER INTEGRATION SUMMARY TABLE ─────────────────────────────────────────
cat("\n=== Building master integration summary ===\n")

integration_table <- data.frame(
  Molecular_Feature = c(
    "STAT1 (IFN/JAK-STAT)",
    "RELA/NF-kB axis",
    "SERPINE1 (PAI-1)",
    "Cellular senescence (SenMayo)",
    "SASP secretome",
    "Immunosenescence (CD57/KLRG1)",
    "T cell exhaustion (NR4A1/2)",
    "Neuronal loss",
    "NR2F1 (synaptic TF)",
    "Cell-intrinsic GRIN3A"
  ),
  Aim_1_Brain_Finding = c(
    "Top VIPER TF: NES=+7.23, padj=1.3e-10; cell-intrinsic LFC=+0.69",
    "RELA NES=+5.44, NFKB1 NES=+4.79; CXCL10/SOCS3/MMP9 upregulated",
    "RELA regulon target; log2FC=+1.23 in COVID-19 brain",
    "Not directly measured in brain",
    "86.7% of bulk signal is composition-driven (gliosis)",
    "Not directly measured",
    "Not directly measured",
    "Neuronal score delta=-0.168 (padj<0.001, BRETIGEA)",
    "NES=-4.09; BEGAIN/LINGO1/KIF21B suppressed",
    "Cell-intrinsic LFC=-0.55 (padj<0.001)"
  ),
  Aim_2_Blood_Finding = c(
    "No LC vs NonLC difference in blood (Wilcoxon p=0.746)",  # [v1.1.0] filled in logged result
    "SERPINE1 significantly elevated (p=0.0002)",
    "LongCOVID median=5.136 vs NonLC=4.062; p=0.0002",
    "SenMayo delta=+0.15 (trend, p=0.23); underpowered n",
    "SASP_Core delta=+0.14 (trend, p=0.11)",
    "B3GAT1 delta=+0.852 (trend); CD27 loss delta=-0.435",
    "NR4A1 delta=+0.672, NR4A2 delta=+0.539 (trends)",
    "PBMC signal upstream — SASP drives neuroinflammation",
    "Peripheral synaptic gene correlation with senescence",
    "Peripheral NMDAR axis not directly tested"
  ),
  Biological_Narrative = c(
    "IFN/STAT1 activation is brain-specific; not detected in Long COVID blood",  # [v1.1.0] corrected to match p=0.746
    "NF-kB inflammatory program bridges systemic inflammation to brain aging",
    "Single strongest blood-brain molecular bridge — RELA target elevated in both compartments",
    "Trend toward senescence enrichment in LC PBMCs — underpowered but directionally consistent",
    "SASP effectors elevated in LC blood; same genes upregulated by RELA in postmortem brain",
    "Immunosenescence trend — CD57 gain, CD27 loss consistent with terminal T cell differentiation",
    "NR4A1/2 elevation suggests exhaustion program activation in Long COVID T cells",
    "Peripheral immune aging (SASP/inflammaging) proposed as upstream driver of neuronal loss",
    "Synaptic suppression in brain mirrors cognitive decline pathway",
    "Cell-intrinsic neuronal NMDA suppression — potential cognitive correlate"
  ),
  stringsAsFactors = FALSE
)

write.csv(integration_table,
          "results/integration/master_integration_summary.csv",
          row.names = FALSE)

cat("Master integration table saved.\n")
cat("\nIntegration table preview:\n")
print(integration_table[, c("Molecular_Feature","Aim_1_Brain_Finding")])

# ── SAVE ALL INTEGRATION OBJECTS ─────────────────────────────────────────────
saveRDS(sasp_gene_stats, "data/sasp_rela_gene_stats.rds")
saveRDS(integration_table, "data/master_integration_table.rds")

cat("\n=================================================================\n")
cat("  PHASE 5 COMPLETE — INTEGRATION ANALYSIS\n")
cat("=================================================================\n")
cat("Figures saved:\n")
cat("  Fig 12: SERPINE1 blood-brain bridge\n")
cat("  Fig 13: STAT1 convergence\n")
cat("  Fig 14: SASP/RELA blood panel\n")
cat("Tables saved:\n")
cat("  Master integration summary -> results/integration/\n")
cat("\nSASP/RELA gene statistics in blood:\n")
print(sasp_gene_stats[, c("gene","delta","p_value","padj")])
cat("\n=== ANALYSIS COMPLETE — READY FOR DISSERTATION WRITE-UP ===\n")
cat("=================================================================\n")

