# SERPINE1/PAI-1 as a Cross-Compartment Biomarker of Long COVID–Associated Cognitive Aging

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23754758.svg)](https://doi.org/10.5281/zenodo.23754758)

R analysis code for the manuscript (in preparation):

> **Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with NF-κB Activation in COVID-19 Postmortem Brain: A Cross-Compartment Transcriptomic Biomarker of Cognitive Aging**
> Kelechi Wisdom Elechi · Gabriel A. de Erausquin
> Glenn Biggs Institute for Alzheimer's & Neurodegenerative Diseases · UT Health San Antonio · Barshop Institute

---

## Version note (v1.1.0)

v1.0.0 contained scripts rewritten after the fact; they did not match the analysis that produced the reported results (different sample grouping, normalization, gene sets and correlation method). **v1.1.0 replaces them with the code that was actually executed**, taken from the R console logs of that run. Every change from the logged code is marked `# [v1.1.0]` in the scripts and listed in each script header. Use v1.1.0 or later.

---

## Datasets

| GEO accession | Description | Used as |
|---|---|---|
| [GSE224615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224615) | Long COVID vs. non-Long COVID PBMC bulk RNA-seq (Yin et al., 2023, *Nature Immunology*) | Primary analysis (36 samples; 35 after outlier removal: LC = 22, non-LC = 13) |
| [GSE188847](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188847) | COVID-19 vs. control postmortem frontal cortex RNA-seq | Brain values taken from the Aim 1 analysis (VIPER/DoRothEA, DESeq2) and hard-coded in script 04 with their source cited |

**Before running:** download the supplementary file `GSE224615_DEGs.xlsx` from the GSE224615 GEO page and place it in `data/`. Script 01 downloads the series-matrix metadata itself with GEOquery.

---

## Pipeline

Run from the repository root, in order:

```bash
Rscript scripts/01_data_qc.R                  # metadata, expression matrix, QC, PCA, outlier flags
Rscript scripts/02_senescence_scoring.R       # removes GSM7027490; GSVA senescence scores
Rscript scripts/03_immunosenescence.R         # immunosenescence GSVA, marker genes, co-enrichment
Rscript scripts/04_integration.R              # SERPINE1 bridge, STAT1, SASP/RELA panel, master table
Rscript scripts/05_effect_sizes_sensitivity.R # effect sizes + sensitivity analyses (new in v1.1.0)
```

| Script | What it does |
|---|---|
| `01_data_qc.R` | Reads GEO metadata (`lc status:ch1`, sex, pid); builds the 18,040-gene × 36-sample matrix from `GSE224615_DEGs.xlsx` (`Gene.name` symbols, duplicates removed); `log2(x + 1)` transform; zero-variance filter; PCA on all genes; flags samples > 3 SD from their group centroid on PC1/PC2 |
| `02_senescence_scoring.R` | Removes outlier GSM7027490 (LC; 4.42 SD on PC1; GSM7027478 at 3.31 SD retained); GSVA (`gsvaParam`, default Gaussian kernel, minSize 5, maxSize 500) for SenMayo (57 genes), CellAge (46), SASP Core (12), p21/p16/p53 (7); two-sided Wilcoxon rank-sum (group medians reported), BH across the 4 sets |
| `03_immunosenescence.R` | GSVA for Immunosenescence Panel (13), Exhaustion TFs (10), Exhaustion Signature (29), NK Immunosenescence (8); Wilcoxon on 15 marker genes (BH); Pearson co-enrichment between score sets |
| `04_integration.R` | Wilcoxon for a 9-gene SASP/RELA panel (BH across the 9 genes); SERPINE1 blood–brain bridge figure; STAT1 blood test; master integration table |
| `05_effect_sizes_sensitivity.R` | Rank-biserial r for every Wilcoxon test; Hodges–Lehmann shift for SERPINE1; SERPINE1 in the transcriptome-wide DESeq2 table; SERPINE1 with the outlier included; Spearman version of the co-enrichment test; Supplementary Figure S1 (Immunosenescence Panel vs SASP Core scatter) |

---

## Key results (from the logged run, n = 35)

| Analysis | Result |
|---|---|
| SERPINE1, LC vs non-LC (9-gene SASP/RELA panel) | median 5.136 vs 4.062; Δ = +1.073; Wilcoxon p = 0.00021; BH padj = 0.0019 (9 tests)* |
| TNF (same panel) | Δ = +0.489; p = 0.012; padj = 0.054 |
| Senescence GSVA (median Δ, p) | SenMayo +0.152, 0.226 · CellAge +0.090, 0.506 · SASP Core +0.143, 0.112 · p21/p16/p53 −0.004, 0.573 |
| Immunosenescence GSVA | no set significant (all padj = 0.986) |
| Immunosenescence Panel vs SASP Core | Pearson r = −0.337, p = 0.047 |
| STAT1 in blood | Wilcoxon p = 0.746 |
| SERPINE1 in brain (Aim 1, GSE188847) | log2FC = +1.23, RELA regulon target (RELA NES = +5.44) |

\* Script 04 prints padj = 0.0018 because it rounds p to 4 decimals before the BH step; the exact value, from unrounded p, is 0.0019 (reported by script 05).

Only 2 genes (OR7D2, ALAS2) reach padj < 0.05 in the transcriptome-wide pre-computed DESeq2 comparison; the SERPINE1 result is from a hypothesis-driven 9-gene panel. Script 05 reports SERPINE1's transcriptome-wide values alongside it.

---

## Software

R 4.3.2 (macOS). Bioconductor: GEOquery 2.70, GSVA 1.50 (versions as recorded at analysis time; `sessionInfo()` output to be added). CRAN: readxl, dplyr, tidyr, ggplot2 3.5, ggrepel, ggbeeswarm, patchwork, pheatmap, RColorBrewer, viridis, scales.

```r
install.packages(c("readxl","dplyr","tidyr","ggplot2","ggrepel","ggbeeswarm",
                   "patchwork","pheatmap","RColorBrewer","viridis","scales"))
BiocManager::install(c("GEOquery","GSVA","org.Hs.eg.db"))
```

---

## Citation

> Elechi KW, de Erausquin GA (2026). *SERPINE1 Long COVID Aging Analysis Pipeline*. Zenodo. https://doi.org/10.5281/zenodo.23754758

## Contact

Kelechi Wisdom Elechi — kelechiwisdom20@gmail.com

## License

Code: MIT. Data: subject to the GEO terms for GSE224615 and GSE188847.
