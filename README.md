# SERPINE1/PAI-1 as a Cross-Compartment Biomarker of Long COVID–Associated Cognitive Aging

Reproducible R analysis pipeline for the manuscript:

> **Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with NF-κB Activation in COVID-19 Postmortem Brain: A Cross-Compartment Transcriptomic Biomarker of Cognitive Aging**  
> Kelechi Wisdom Elechi · Gabriel A. de Erausquin  
> Glenn Biggs Institute for Alzheimer's & Neurodegenerative Diseases · UT Health San Antonio · Barshop Institute  
> *Submitted to Frontiers in Aging Neuroscience*

---

## Overview

This repository contains all analysis code for a cross-compartment transcriptomic study linking peripheral SASP/senescence gene expression in Long COVID PBMCs (blood, **Aim 2**) to NF-κB–driven neuroinflammation in COVID-19 postmortem brain (**Aim 1**).

**Primary finding:** SERPINE1 (PAI-1) is the sole SASP Core gene passing FDR correction in Long COVID peripheral blood mononuclear cells (Wilcoxon p = 0.0002, adj p = 0.0018), and is independently upregulated in COVID-19 postmortem frontal cortex (log2FC = +1.23, padj = 0.024), establishing a cross-compartment convergence signal.

---

## Datasets

| GEO Accession | Description | n (post-QC) |
|---|---|---|
| [GSE224615](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224615) | Long COVID PBMC bulk RNA-seq (Yin et al., 2023, *Nature Immunology*) | 35 (LC=22, NonLC=13) |
| [GSE188847](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188847) | COVID-19 postmortem frontal cortex RNA-seq | 51 |

**Data access:** Raw count files must be downloaded manually from GEO and placed in `data/`:
- `data/GSE224615_DEGs.xlsx` — from GSE224615 series matrix
- `data/GSE188847_deseq2_results.csv` — DESeq2 output from Aim 1 pipeline
- `data/GSE188847_viper_nes.csv` — VIPER/DoRothEA NES table from Aim 1 pipeline

---

## Repository Structure

```
├── scripts/
│   ├── 01_download_qc.R          # Data loading, QC, normalization, PCA
│   ├── 02_sasp_spotlight.R       # SASP gene panel + SERPINE1 primary analysis
│   ├── 03_senescence_scoring.R   # GSVA senescence scoring (5 gene sets)
│   ├── 04_immunosenescence.R     # Immunosenescence profiling + co-enrichment
│   └── 05_integration_bridge.R  # Cross-compartment integration + bridge figures
├── data/                         # Raw data (not tracked — see Data access above)
├── results/
│   ├── figures/                  # Output figures (not tracked)
│   └── *.rds / *.csv             # Intermediate results (not tracked)
└── README.md
```

---

## Analysis Pipeline

Run scripts in order from the project root directory:

```r
source("scripts/01_download_qc.R")          # ~1 min
source("scripts/02_sasp_spotlight.R")        # ~2 min
source("scripts/03_senescence_scoring.R")    # ~5 min (GSVA)
source("scripts/04_immunosenescence.R")      # ~2 min
source("scripts/05_integration_bridge.R")   # ~2 min
```

Or from the shell:

```bash
Rscript scripts/01_download_qc.R
Rscript scripts/02_sasp_spotlight.R
Rscript scripts/03_senescence_scoring.R
Rscript scripts/04_immunosenescence.R
Rscript scripts/05_integration_bridge.R
```

Each script is self-contained and documents its inputs/outputs in the header.

---

## Key Results

| Analysis | Statistic | Significance |
|---|---|---|
| SERPINE1 blood (Wilcoxon) | LC median=5.136 vs NonLC=4.062 | p=0.0002, **adj p=0.0018** ✓ |
| SERPINE1 brain (DESeq2) | log2FC=+1.23 | padj=0.024 ✓ |
| SenMayo GSVA | delta=+0.152 | p=0.226 (trend) |
| SASP Core GSVA | delta=+0.143 | p=0.112 (trend) |
| ImmunoSen–SASP co-enrichment | Spearman r=−0.337 | p=0.047 ✓ |
| STAT1 blood | — | p=0.746 (null; constitutive) |
| RELA brain VIPER | NES=+5.44 | p<1×10⁻⁸ ✓ |
| NFKB1 brain VIPER | NES=+4.79 | p<1×10⁻⁷ ✓ |

---

## Software Requirements

```r
# CRAN
install.packages(c("readxl", "dplyr", "tidyr", "ggplot2", "ggrepel",
                   "ggpubr", "patchwork", "pheatmap", "RColorBrewer",
                   "readr"))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("GSVA")

# Optional (for Aim 1 brain analysis)
BiocManager::install(c("viper", "dorothea", "DESeq2", "BRETIGEA"))
```

Tested with R ≥ 4.3.0. GSVA ≥ 1.50.0 uses the `gsvaParam()` API (Bioc 3.18+).

---

## Figure Map

| Figure | Script | Description |
|---|---|---|
| Fig 1 | `01_download_qc.R` | PCA of GSE224615 samples |
| Fig 2 | `02_sasp_spotlight.R` | SASP / NF-κB 9-gene panel |
| Fig 3 | `03_senescence_scoring.R` | GSVA scores — 5 senescence gene sets |
| Fig 4 | `04_immunosenescence.R` | ImmunoSen–SASP co-enrichment scatter |
| Fig 5 | `05_integration_bridge.R` | SERPINE1 blood-brain bridge |
| Fig 6 | `05_integration_bridge.R` | STAT1 compartment divergence |

---

## Citation

If you use this code, please cite:

> Elechi KW, de Erausquin GA. Peripheral SERPINE1/PAI-1 Elevation in Long COVID Converges with NF-κB Activation in COVID-19 Postmortem Brain: A Cross-Compartment Transcriptomic Biomarker of Cognitive Aging. *Frontiers in Aging Neuroscience* (submitted 2026).

---

## Contact

Kelechi Wisdom Elechi — kelechiwisdom20@gmail.com  
Glenn Biggs Institute for Alzheimer's & Neurodegenerative Diseases  
UT Health San Antonio · Barshop Institute for Longevity and Aging Studies

---

## License

Code: MIT License  
Data: Subject to GEO data use agreements for GSE224615 and GSE188847.
