# Shotgun metagenomic analysis reveals taxonomic and functional alterations in the gut microbiome across prodromal and symptomatic Lewy body disease

**DOI:** 

**Authors:** Xiaowei Zhao, Stuart J. McCarter, Vinod K. Gupta, Kiera M. Grant, Erik K. St. Louis, Kejal Kantarci, Rodolfo Savica, Max Hill, Helen E. Vuong, Christopher Staley, Bradley F. Boeve, Owen A. Ross, Levi M. Teigen, and Jaeyun Sung

## Overview

![study_design](./study_design.jpg)


This repository contains the R analysis code for a gut microbiome study comparing individuals with Lewy Body Disease (LBD) and isolated REM Sleep Behavior Disorder (iRBD) against their matched household controls. Shotgun metagenomic data were profiled using MetaPhlAn (taxonomy) and HUMAnN (functional pathways and gene families).

## Study Design

Four groups are compared across analyses:

| Group | Label |
|---|---|
| Lewy Body Disease patients | `LBD` |
| LBD cohabitant controls | `LBD-Control` |
| Isolated REM Sleep Behavior Disorder patients | `iRBD` |
| iRBD household controls | `iRBD-Control` |

Comparisons: LBD vs. LBD-Control and iRBD vs. iRBD-Control

## Repository Structure

```
├── 0-raw_data/
│   ├── metaphlan_results_new.tsv      # Raw relative abundance data (MetaPhlAn output)
│   ├── pathway_results.tsv            # Functional pathway abundance (HUMAnN output)
│   ├── imputed_BMI_metadata.csv       # Sample metadata (condition, age, sex, BMI, clinical scores)
│   └── taxonomy_info                  # Taxonomy information at each level
│
├── 1-alpha_beta_diversity/
│   ├── alpha_beta_diversity_final_version.R
│   ├── shannon_4_groups.pdf                        # Shannon index boxplots across all 4 groups
│   ├── richness_4_groups.pdf                       # Species richness boxplots across all 4 groups
│   ├── stacked_bar_plot_family_all.pdf             # Stacked bar plots of family-level relative abundance
│   ├── distribution_family_relative_abundance.csv  # Family-level relative abundance table
│   ├── lbd_vs_control_results/
│   │   ├── alpha_diversity_res_lbd_vs_control.csv  # Alpha diversity metrics for LBD vs. LBD-Control
│   │   └── pcoa_lbd_vs_control.pdf                 # PCoA plot for LBD vs. LBD-Control
│   └── irbd_vs_control_results/
│       ├── alpha_diversity_results_irbd_vs_control.csv  # Alpha diversity metrics for iRBD vs. iRBD-Control
│       └── pcoa_irbd_vs_control.pdf                     # PCoA plot for iRBD vs. iRBD-Control
│
├── 2-microbial_taxa_analysis/
│   ├── species_preprocessed.csv                            # Preprocessed species relative abundance (remove low relateve abundance data, output of script 1, input for script 2)
│   ├── 1-differential_abundance_analysis_final_version.R
│   ├── 1-lbd_vs_control_diff_abundance_results/
│   │   ├── lbd_vs_control_differential_abundance_p_values.csv   # Differential abundance analysis results with P- and q-values (LBD vs. LBD-Control)
│   │   ├── lbd_vs_control_diff_abun_p_fc_done.csv               # Differential abundance analysis results with P- and q-values, as well as log2FC (LBD vs. LBD-Control)
│   │   ├── sig_lbd_vs_control_diff_abun_species.csv             # Significant differential abundance species
│   │   ├── species_higher_in_LBD_than_Control.csv
│   │   ├── species_higher_in_Control_than_LBD.csv
│   │   ├── cohensD_lbd_vs_control_species.csv                   # Cohen's D effect sizes
│   │   ├── volcano_plot_LBD_vs_Control.pdf
│   │   ├── lbd_vs_control_boxplot_species.pdf
│   │   └── cohensD_LBD_species.pdf
│   ├── 2-different_prevalence_analysis_final_version.R
│   ├── 2-irbd_vs_control_diff_abundance_results/
│   │   ├── irbd_vs_control_differential_abundance_p_values.csv  # Differential abundance analysis results with P- and q-values (iRBD vs. iRBD-Control)
│   │   ├── irbd_vs_control_diff_abun_p_fc_done.csv              # Differential abundance analysis results with P- and q-values, as well as log2FC (iRBD vs. iRBD-Control)
│   │   ├── species_higher_in_iRBD_than_Control.csv
│   │   ├── species_higher_in_Control_than_iRBD.csv
│   │   ├── cohensD_irbd_vs_control_species.csv                  # Cohen's D effect sizes
│   │   ├── volcano_plot_iRBD_vs_Control.pdf
│   │   ├── irbd_vs_control_boxplot_species.pdf
│   │   └── cohensD_irbd_vs_control_species.pdf
│   ├── 3-lbd_vs_control_diff_prevalence_results/
│   │   ├── lbd_vs_control_prevalence_results.csv
│   │   └── sig_lbd_vs_control_prev.csv                          # Significant differentially prevalent species
│   └── 4-irbd_vs_control_diff_prevalence_results/
│       ├── irbd_vs_control_prevalence_results.csv
│       └── sig_irbd_vs_control_prev.csv                         # Significant differentially prevalent species
│
├── 3-microbial_functional_pathway_analysis/
│   ├── pathway_analysis_final_version.R
│   ├── pathway_taxonomy_analysis_final_version.R
│   └── pathways_contributed_by_bacterial_species_final_version.R
│
├── 4-gene_families_analysis/
│   └── gene_families_analysis_final_version.R
│
└── 5-correlation_analysis/
    └── partial_correlation_analysis.R
```

## Analysis Modules

### 1. Alpha & Beta Diversity (`1-alpha_beta_diversity/`)

- Computes Shannon index, Simpson index, inverse Simpson index, and species richness per sample
- Statistical testing via mixed-effects linear models (`lmerTest`) with household as a random effect, and Wilcoxon rank-sum tests
- Beta diversity via Bray-Curtis dissimilarity with arcsine square-root transformation; PERMANOVA (`vegan::adonis2`) stratified by household; PCoA plots with 95% confidence ellipses (`ade4`, `ggplot2`)
- Stacked bar plots of family-level relative abundance across all four groups

### 2. Microbial Taxa Analysis (`2-microbial_taxa_analysis/`)

**Script 1 — Differential Abundance** (`1-differential_abundance_analysis_final_version.R`):
- Species with relative abundance below 10^−4.7 are set to 0
- Species detected in fewer than 10% of samples excluded
- Relative abundances arcsine square-root transformed prior to statistical testing
- Differentially abundant species were identified using mixed-effects linear models (`lmerTest`) with household as a random effect and age and BMI as covariates
- Volcano plots, per-species boxplots, and Cohen's D effect sizes for both LBD vs. LBD-Control and iRBD vs. iRBD-Control

**Script 2 — Differential Prevalence** (`2-different_prevalence_analysis_final_version.R`):
- Species with relative abundance below 10^−4.7 are set to 0
- Species detected in fewer than 10% of samples excluded
- Relative abundances arcsine square-root transformed prior to statistical testing
- Differential prevalence analysis using Fisher's exact test
- Results saved separately for LBD vs. LBD-Control and iRBD vs. iRBD-Control

### 3. Functional Pathway Analysis (`3-microbial_functional_pathway_analysis/`)

- Mixed-effects linear models to identify differentially abundant metabolic pathways
- Volcano plots of pathway effect sizes vs. significance
- Attribution of pathway abundance to contributing bacterial species (HUMAnN stratified output)

### 4. Gene Families Analysis (`4-gene_families_analysis/`)

- Differential abundance analysis of UniRef90 gene families (HUMAnN output)
- Abundance cutoff (10^−7) applied to filter low-abundance gene families
- Mixed-effects linear models with household as a random effect

### 5. Partial Correlation Analysis (`5-correlation_analysis/`)

- Spearman partial correlations between microbial features (species and pathways) and clinical measurements (e.g., MoCA, CDR, UPDRS)
- Covariates: age, BMI
- Two-stage screen: fast parametric `pcor.test` followed by permutation-based validation

## R Package Requirements

```r
install.packages(c(
  "lmerTest", "vegan", "ade4", "ggplot2", "reshape2",
  "dplyr", "ggpubr", "tidyr", "RColorBrewer",
  "pheatmap", "ppcor", "tidyverse", "ggrepel", "purrr",
  "effectsize", "readxl", "writexl"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
```

## Data Preprocessing

All scripts apply a consistent preprocessing pipeline before statistical testing:

1. Extract species-level rows from MetaPhlAn output (rows matching `s__` but not `t__`)
2. Normalize counts to relative proportions (column sums = 1)
3. Apply an abundance cutoff (features below this threshold in a sample are set to zero)
4. Remove features detected in fewer than 10% of samples 
5. Arcsine square-root transformation for beta diversity and differential abunance analyses
6. Merge with sample metadata; analyses use matched case-control pairs (household ID as random effect)

## Reproducibility

PERMANOVA analyses use `set.seed(10)` before each permutation test to ensure reproducible results.
