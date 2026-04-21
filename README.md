# Cui MADA Project — HMP 16S Microbiome Analysis Across Five Major Body Sites

This repository contains the full data-analysis project for the MADA course.
The project uses 16S rRNA V3–V5 data from the Human Microbiome Project (HMP)
and tests whether alpha diversity, community composition, and microbial
co-occurrence network structure differ across five major body sites:
Airways, Gut, Oral, Skin and Urogenital.

Main research questions:

- **H1** — Alpha diversity differs across the five body sites.
- **H2** — Samples cluster by body site in beta-diversity ordination space.
- **H3** — Co-occurrence networks differ in structure across body sites.

## Repository structure

```
Cui-MADA-project/
├── assets/                       # Bibliography, CSL, schematics
├── code/
│   ├── processing-code/          # HMP data download + cleaning (step 1)
│   ├── eda-code/                 # Exploratory data analysis (step 2)
│   └── analysis-code/            # Formal stats + networks (step 3)
├── data/
│   ├── raw-data/                 # Inputs (HMP is downloaded at runtime)
│   └── processed-data/           # ps_filt.rds, ps_rel.rds
├── products/
│   ├── manuscript/               # Manuscript.qmd + supplement
│   ├── report/                   # HTML report
│   ├── presentation/             # Slides
│   └── poster/                   # Poster
└── results/
    ├── figures/                  # All figures used in the manuscript
    ├── tables/                   # All tables used in the manuscript (RDS)
    └── output/                   # Larger intermediate objects (networks)
```

## Prerequisites

- R ≥ 4.3 and Quarto
- The following R packages (CRAN + Bioconductor):

  ```r
  install.packages(c(
    "tidyverse", "here", "broom", "vegan", "igraph", "ggraph",
    "tidygraph", "Hmisc", "RColorBrewer", "skimr", "readxl",
    "knitr", "rmarkdown"
  ))

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c(
    "phyloseq", "HMP16SData", "SummarizedExperiment", "microbiome"
  ))
  ```

## Reproduction instructions

Run the following three steps in order. Each step writes its outputs to
`data/processed-data/` or `results/`. The working directory is set
automatically by the `here` package so commands can be run from anywhere in
the project.

1. **Processing** — downloads the HMP V3–V5 data via `HMP16SData::V35()`,
   converts it to a `phyloseq` object, removes zero-count taxa, filters
   rare taxa (present in ≤ 5 samples), transforms to relative abundance,
   and writes `data/processed-data/ps_filt.rds` and
   `data/processed-data/ps_rel.rds`:

   ```r
   quarto::quarto_render(
     here::here("code", "processing-code", "processingfile-v1.qmd"))
   ```

2. **Exploratory data analysis** — renders sequencing-depth, top-phyla,
   alpha-diversity and PCoA-preview figures, and writes the
   `exploratory_summary_by_subsite.rds` table:

   ```r
   quarto::quarto_render(here::here("code", "eda-code", "eda.qmd"))
   ```

3. **Formal statistical analysis** — runs Kruskal-Wallis + pairwise
   Wilcoxon for alpha diversity, PERMANOVA + BETADISPER for beta diversity,
   and builds per-site co-occurrence networks. Writes all tables and
   figures used in the manuscript:

   ```r
   source(here::here("code", "analysis-code", "statistical-analysis.R"))
   ```

4. **Manuscript and supplement** — rendering these files pulls in the
   saved tables and figures:

   ```r
   quarto::quarto_render(
     here::here("products", "manuscript", "Manuscript.qmd"))
   quarto::quarto_render(
     here::here("products", "manuscript", "supplement",
                "Supplementary-Material.qmd"))
   ```

## File conventions

- Lower-case names with `-` as the word separator.
- R scripts end in `.R`; Quarto files in `.qmd`.
- Tables are saved as `.rds`; figures as `.png`.
- All randomness uses `set.seed(123)` for reproducibility.

## Data

The primary input is the HMP 16S V3–V5 dataset obtained at runtime via
`HMP16SData::V35()`. No raw files need to be stored in the repository; the
download is reproducible and checked against the package version. The small
`exampledata.xlsx` under `data/raw-data/` is left over from the template and
is not used in this analysis.

## Citations

Example reproducible projects referenced here include McKay et al.
(2020a, 2020b); see `assets/dataanalysis-references.bib`.
