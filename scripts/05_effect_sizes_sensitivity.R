# =============================================================================
# Script 05: Effect sizes and sensitivity analyses  [new in v1.1.0]
# Project: Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with
#          NF-kB Activation in COVID-19 Postmortem Brain
# Author:  Kelechi Wisdom Elechi
# -----------------------------------------------------------------------------
# Adds reporting only; it does not change any result from scripts 01-04.
#   1. Rank-biserial correlation (r_rb) for every Wilcoxon comparison
#      r_rb = 2W / (n_LC * n_NonLC) - 1 ; ranges -1..1, positive = higher in LC
#   2. Hodges-Lehmann location shift with 95% CI for SERPINE1
#   3. SERPINE1 in the transcriptome-wide pre-computed DESeq2 results
#   4. Sensitivity: SERPINE1 / SASP-RELA panel with outlier GSM7027490 included
#   5. Sensitivity: co-enrichment correlation with Spearman's rho
#   6. Supplementary Figure S1: Immunosenescence Panel vs SASP Core scatter
# INPUT:  data/*.rds from scripts 01-04
# OUTPUT: results/integration/effect_sizes_all_tests.csv
#         results/integration/sensitivity_outlier_included.csv
#         results/integration/sensitivity_summary.txt
#         results/integration/figS1_coenrichment_sasp_core.png/.pdf
# =============================================================================

suppressPackageStartupMessages({ library(dplyr) })
dir.create("results/integration", showWarnings = FALSE, recursive = TRUE)

expr_clean  <- readRDS("data/expr_outlier_removed.rds")      # 35 samples
col_clean   <- readRDS("data/coldata_outlier_removed.rds")
expr_all    <- readRDS("data/expr_log2_symbol_gse224615.rds") # 36 samples
col_all     <- readRDS("data/coldata_gse224615.rds")
sen_scores  <- readRDS("data/gsva_senescence_scores.rds")
isen_scores <- readRDS("data/gsva_immunosenescence_scores.rds")
coenrich_df <- readRDS("data/coenrichment_scores.rds")
deg_results <- readRDS("data/deseq2_results_LC_vs_NonLC.rds")

wilcox_es <- function(lc, nlc, label, family) {
  lc <- as.numeric(lc); nlc <- as.numeric(nlc)
  wt <- wilcox.test(lc, nlc, exact = FALSE)
  data.frame(
    family       = family,
    feature      = label,
    n_LC         = length(lc),
    n_NonLC      = length(nlc),
    median_LC    = median(lc),
    median_NonLC = median(nlc),
    delta        = median(lc) - median(nlc),
    W            = unname(wt$statistic),
    p_value      = wt$p.value,
    r_rb         = unname(2 * wt$statistic / (length(lc) * length(nlc)) - 1),
    stringsAsFactors = FALSE
  )
}

lc_ids  <- col_clean$sample_id[col_clean$group == "LongCOVID"]
nlc_ids <- col_clean$sample_id[col_clean$group == "Non_LongCOVID"]

# ── 1. Effect sizes for every test family (BH within family, as in 02-04) ────
panel_genes  <- intersect(c("SERPINE1","IL6","MMP9","CCL2","CXCL10",
                            "VEGFA","ICAM1","IL1B","CXCL8","TNF"),
                          rownames(expr_clean))
marker_genes <- intersect(c("KLRG1","B3GAT1","TIGIT","LAG3","PDCD1",
                            "TOX","NR4A1","NR4A2","CD28","CD27",
                            "SELL","IL7R","TCF7","FCGR3A","GZMB"),
                          rownames(expr_clean))

by_gene <- function(genes, family)
  bind_rows(lapply(genes, function(g)
    wilcox_es(expr_clean[g, lc_ids], expr_clean[g, nlc_ids], g, family)))
by_set <- function(scores, family)
  bind_rows(lapply(rownames(scores), function(s)
    wilcox_es(scores[s, lc_ids], scores[s, nlc_ids], s, family)))

es <- bind_rows(
  by_gene(panel_genes,  "SASP/RELA 9-gene panel"),
  by_set(sen_scores,    "Senescence GSVA"),
  by_set(isen_scores,   "Immunosenescence GSVA"),
  by_gene(marker_genes, "Immunosenescence markers"),
  wilcox_es(expr_clean["STAT1", lc_ids], expr_clean["STAT1", nlc_ids],
            "STAT1", "STAT1 (single test)")
) %>%
  group_by(family) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

write.csv(es, "results/integration/effect_sizes_all_tests.csv", row.names = FALSE)
cat("\n=== Effect sizes (rank-biserial r) ===\n")
print(as.data.frame(es %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))

# ── 2. Hodges-Lehmann shift for SERPINE1 ─────────────────────────────────────
hl <- wilcox.test(as.numeric(expr_clean["SERPINE1", lc_ids]),
                  as.numeric(expr_clean["SERPINE1", nlc_ids]),
                  exact = FALSE, conf.int = TRUE)

# ── 3. SERPINE1 in the transcriptome-wide DESeq2 results ─────────────────────
gene_col <- intersect(c("gene", "Gene.name", "symbol"), colnames(deg_results))[1]
serp_deg <- deg_results[deg_results[[gene_col]] == "SERPINE1", ]

# ── 4. Sensitivity: include outlier GSM7027490 (36 samples) ──────────────────
lc_all  <- col_all$sample_id[col_all$group == "LongCOVID"]
nlc_all <- col_all$sample_id[col_all$group == "Non_LongCOVID"]
sens <- bind_rows(lapply(panel_genes, function(g)
  wilcox_es(expr_all[g, lc_all], expr_all[g, nlc_all], g,
            "SASP/RELA panel, outlier included"))) %>%
  mutate(padj = p.adjust(p_value, method = "BH"))
write.csv(sens, "results/integration/sensitivity_outlier_included.csv", row.names = FALSE)

# ── 5. Sensitivity: Spearman co-enrichment ───────────────────────────────────
sp <- cor.test(coenrich_df$ImmunoSen, coenrich_df$SASP_Core,
               method = "spearman", exact = FALSE)
pe <- cor.test(coenrich_df$ImmunoSen, coenrich_df$SASP_Core, method = "pearson")

# ── 6. Supplementary Figure S1: Immunosenescence Panel vs SASP Core ─────────
# (Script 03's fig11 plots ImmunoSen vs SenMayo; this plots the SASP Core pair
#  whose correlation, r = -0.337, is reported in the text.)
suppressPackageStartupMessages({ library(ggplot2) })
grp_cols <- c("Non_LongCOVID" = "#2196F3", "LongCOVID" = "#F44336")
cd <- coenrich_df
if (!"group" %in% colnames(cd))
  cd$group <- col_clean$group[match(cd$sample_id, col_clean$sample_id)]
p_s1 <- ggplot(cd, aes(x = ImmunoSen, y = SASP_Core, colour = group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE,
              colour = "black", linetype = "dashed", linewidth = 0.8) +
  scale_colour_manual(values = grp_cols, name = "Group") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = sprintf("Pearson r = %.3f, p = %.3f\nSpearman rho = %.3f, p = %.3f",
                           pe$estimate, pe$p.value, sp$estimate, sp$p.value)) +
  labs(title = "Co-Enrichment: Immunosenescence Panel vs SASP Core",
       subtitle = "GSE224615 | all 35 samples",
       x = "Immunosenescence Panel Score (GSVA)",
       y = "SASP Core Score (GSVA)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
ggsave("results/integration/figS1_coenrichment_sasp_core.png", p_s1,
       width = 7, height = 5.5, dpi = 300)
ggsave("results/integration/figS1_coenrichment_sasp_core.pdf", p_s1,
       width = 7, height = 5.5)

# ── Summary for the manuscript ───────────────────────────────────────────────
s1  <- es[es$feature == "SERPINE1", ]
s1o <- sens[sens$feature == "SERPINE1", ]
sink("results/integration/sensitivity_summary.txt")
cat("SERPINE1 (35 samples, primary analysis)\n")
cat(sprintf("  median LC = %.3f, NonLC = %.3f, delta = %.3f\n",
            s1$median_LC, s1$median_NonLC, s1$delta))
cat(sprintf("  Wilcoxon W = %.0f, p = %.2g, BH padj (9-gene panel) = %.2g\n",
            s1$W, s1$p_value, s1$padj))
cat(sprintf("  rank-biserial r = %.2f\n", s1$r_rb))
cat(sprintf("  Hodges-Lehmann shift = %.3f (95%% CI %.3f to %.3f)\n",
            hl$estimate, hl$conf.int[1], hl$conf.int[2]))
cat("\nSERPINE1 in pre-computed transcriptome-wide DESeq2 results (GSE224615_DEGs.xlsx)\n")
print(serp_deg)
cat("\nSERPINE1 with outlier GSM7027490 included (36 samples)\n")
cat(sprintf("  delta = %.3f, p = %.2g, BH padj (9-gene panel) = %.2g, r_rb = %.2f\n",
            s1o$delta, s1o$p_value, s1o$padj, s1o$r_rb))
cat("\nImmunoSen Panel vs SASP Core co-enrichment (35 samples)\n")
cat(sprintf("  Pearson r = %.3f, p = %.3g   (primary, as reported)\n", pe$estimate, pe$p.value))
cat(sprintf("  Spearman rho = %.3f, p = %.3g (sensitivity)\n", sp$estimate, sp$p.value))
sink()
cat("\n"); cat(readLines("results/integration/sensitivity_summary.txt"), sep = "\n")
cat("\nScript 05 complete.\n")
