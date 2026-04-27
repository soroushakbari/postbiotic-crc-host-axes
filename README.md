# Postbiotic-related host transcriptomic axes in colorectal cancer

This repository contains reproducible R code and processed analysis outputs accompanying the manuscript:

**Postbiotic-related host transcriptomic axes in colorectal cancer: evidence mapping and transcriptomic validation across clinical, survival, and tumour microenvironment contexts**

The project evaluates host-side transcriptomic programmes related to major postbiotic-associated biological axes in colorectal cancer, integrating a focused literature evidence map with transcriptomic analyses in TCGA colorectal cancer cohorts and external validation in GSE39582.

## Study overview

The analysis focuses on four postbiotic-related host axes:

- Short-chain fatty acid (SCFA) axis
- Tryptophan / AhR axis
- Polyamine axis
- Ferroptosis / redox axis

The workflow includes:

1. Curation of host gene sets representing postbiotic-related biological axes.
2. Scoring of axis-level transcriptomic activity in colorectal cancer cohorts.
3. Association testing with tumour stage and overall survival.
4. External validation in GSE39582.
5. Correlation of axis scores with tumour microenvironment signatures.
6. Association of axis scores with MSI and immune-high phenotypes.
7. Generation of manuscript-ready tables and figures.

## Repository structure


code/
  R/
    Core R scripts used for gene-set curation, transcriptomic scoring,
    clinical association analyses, survival modelling, TME association
    analyses, and manuscript figure/table generation.

data/
  processed/
    Processed gene-set and figure-input files required to reproduce
    the main analyses and figures.

results/
  tables/
    Final analysis tables and supplementary tables.

  figures_main/
    Final manuscript figures generated from the analysis scripts.

Data sources

This study uses publicly available colorectal cancer transcriptomic and clinical resources, including:

TCGA colorectal cancer cohorts
GSE39582 external validation cohort

Raw public data are not redistributed in this repository. The scripts document the analysis workflow, and processed output files needed for reproducing manuscript tables and figures are provided where appropriate.

Main outputs

The repository includes code and processed outputs for:

Postbiotic-related host axis definitions
Stage association analyses
Univariable and multivariable overall survival models
Tumour microenvironment correlation analyses
MSI and immune-high phenotype associations
Final manuscript figures
Final main and supplementary tables
Reproducibility

The project was developed in R. Required package calls are included in the scripts. Users should run the scripts from the project root directory after adjusting local paths as needed.

Recommended workflow:

source("code/R/02_gene_sets_postbiotics_curate.R")
source("code/R/03_TCGA_download_preprocess.R")
source("code/R/04_TCGA_pathway_scores.R")
source("code/R/05_TCGA_clinical_survival.R")
source("code/R/06_TCGA_TME_signatures.R")
source("code/R/08a_GSE39582_build_expr_clin_v02.R")
source("code/R/08b_GSE39582_scores_assoc_v02.R")
source("code/R/09_OS_multivariable_models.R")
source("code/R/10_TCGA_stratified_MSI_immune_axes.R")
source("code/R/30_make_Figure3_clinical_axes_final_v03.R")
source("code/R/31_make_Supplementary_Table_S4_gene_signature_membership.R")

Figure-generation scripts can be run after the relevant processed input files are available.

Evidence map

The focused evidence map summarizes postbiotic-related evidence in colorectal cancer across mechanistic, clinical, survival, therapeutic-response, and tumour microenvironment contexts. The literature evidence window was updated through April 2026.

License

Code in this repository is released under the MIT License. Public datasets used in this study remain subject to their original access terms and licenses.

Contact

For questions regarding this repository or the accompanying manuscript, please contact:

Soroush Akbari Ardabili
Division of Biochemistry, Department of Basic Sciences
Faculty of Veterinary Medicine, Urmia University
