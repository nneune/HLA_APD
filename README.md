# Repository for all code used for current study

This repository contains the scripts used for the data analysis and figure generation in our study *"Using viral diversity to identify HIV-1 variants under HLA-dependent selection in a systematic viral genome-wide screen."* It includes scripts for generating datasets, processing sequences, calculating APD scores, and creating the figures used in the manuscript. The publication can be found [here](https://doi.org/10.1371/journal.ppat.1012385).

The data generated during the current study cannot be shared publicly due to the sensitive nature and privacy concerns (see [SHCS Open Data Statement](https://www.shcs.ch/294-open-data-statement-shcs)). Investigators with a request for selected data can send a proposal to the [Swiss HIV Cohort Study](https://www.shcs.ch/contact). The provision of data will be evaluated by the Scientific Board of the Swiss HIV Cohort Study and the study team and will be subject to Swiss legal and ethical regulations.

## Main Scripts for Data and Figure Generation

- **`11_data_load.R`**: Loads the primary datasets used in the analysis.
- **`12_manuscript_plots.R`**: Generates all the plots and figures for the manuscript.
- **`13_longitudinal_netMHC.R`**: Handles longitudinal analysis and netMHC predictions and the respective figures.

## Generating FASTA Files and Working with Them

- **`00_formatting_threshhold_to_fasta.R`**: Formats raw NGS data into FASTA files for sequence analysis.
- **`01_sequences_for_alignment.R`**: Prepares sequences for alignment.
- **`01_fastaAA_APD.R`**: Generates amino acid sequences and APD calculations.
- **`01_apd_functions.R`**: Contains functions used for calculating APD scores.

## Additional Functions for Dataset Preparation

- **`10_packages_load.R`**: Loads R packages and functions for analysis.
- **`02_binary_region_aa.R`**: Processes binary data (presence or absence) of amino acids in the sequences.
- **`03_hlaSNPbin_AA.R`**: Prepares binary SNP data for HLA regions.
- **`06_power_Fisher.R`**: Performs power analysis and Fisher's exact test.
- **`forestplot_function.R`**: Generates forest plots for the study.

## Other Data Preparation Scripts

- **`00_HLA_prep.R`**: Cleans and prepares HLA data for analysis.
- **`05_eigensoft.R`**: Performs EIGENSOFT analysis for viral principal components (PCs).
