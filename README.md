# Polygenic Risk Score Simulation and Analysis – Code for Farrer 2026 (Neurobiology of Disease)

This repository contains code and data associated with the manuscript:

**"The challenge of discovering new genes associated with Parkinson’s disease"**  
Farrer, M., Neurobiology of Disease, *Current Opinion in Neurobiology*, 2026 (in press).

## Overview

This repository provides:

- Scripts for the simulation, analysis, and plotting of polygenic risk scores (PRS) as reported in the manuscript, for both "all SNPs" and "known SNPs" analyses.
- An interactive Shiny app for GWAS PRS simulation, available (with identical code) at:  
  **https://bcb-sr.rc.ufl.edu/prs-simulator/**
- Example data used for all analyses and the app.

---

## Repository contents

```
+--- data
|    +--- Leonard_2025_319455_file04.xlsx
|    +--- Leonard_2025_file04_S3.xlsx
+--- R
|    +--- app.R
|    +--- figure-and-summary-stats.R
```

- **data/Leonard_2025_319455_file04.xlsx** and **data/Leonard_2025_file04_S3.xlsx**:  
  Parkinson’s disease GWAS summary statistics, as described below.

- **R/app.R**:  
  The Interactive GWAS Polygenic Risk Score (PRS) Simulator as deployed at UFHCC BCB-SR.  
  This allows users to upload GWAS summary stats, set PRS parameters, and visualize simulated distributions interactively.

- **R/figure-and-summary-stats.R**:  
  Standalone R script that performs the PRS genotype simulation, calculates PRS and OR distributions, creates publication-quality plots (including standardized PRS plots), and generates the summary statistics and figure used in the manuscript.

---

## Data Source

The GWAS summary statistics (data/Leonard_2025_319455_file04.xlsx) originate from:

> **Novel Parkinson’s Disease Genetic Risk Factors Within and Across European Populations**  
> The Global Parkinson’s Genetics Program (GP2), Hampton L. Leonard  
> medRxiv 2025.03.14.24319455; doi: https://doi.org/10.1101/2025.03.14.24319455

*This article is a preprint and has not been peer-reviewed.*

---

## How to use

- To run the interactive PRS simulator, open `R/app.R` in RStudio and click "Run App", or visit https://bcb-sr.rc.ufl.edu/prs-simulator/  
- To reproduce the figures and summary stats from the manuscript, run `R/figure-and-summary-stats.R` after ensuring your working directory includes the `data` folder (with .xlsx files above).

## Acknowledgments

This code and analysis were developed in part from original analysis and foundational code provided by Dr. Heather Kates and the UF Health Cancer Center Biostatistics, Computational Biology, and Shared Resources (UFHCC BCB-SR).

---

For more information, see the manuscript and preprint. For questions, contact hkates@ufl.edu
