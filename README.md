# TCGA-BRCA Bioinformatics Project

## Exploring Molecular Alterations and Their Association with Clinical Characteristics and Overall Survival in Breast Cancer Using TCGA Data

## 1. Project Overview

This project uses publicly available breast cancer data from The Cancer Genome Atlas (TCGA), accessed through cBioPortal, to explore selected genetic alterations and their relationship with clinical characteristics and overall survival.

The analysis focuses on four genes: TP53, BRCA1, BRCA2, and PIK3CA.

The project demonstrates a reproducible bioinformatics workflow using R, covering data acquisition, quality control, data processing, statistical analysis, and survival analysis.

## 2. Research Question

Are selected genetic alterations in breast cancer associated with patient clinical characteristics and overall survival?

## 3. Study Objectives

The objectives of this project are to:

1. Retrieve and process clinical and molecular data from the TCGA breast cancer study.
2. Determine the frequency of mutations in TP53, BRCA1, BRCA2, and PIK3CA.
3. Examine the association between mutation status and tumour stage.
4. Assess differences in overall survival between mutation-status groups.
5. Evaluate the relationship between mutation status and overall survival while accounting for age and tumour stage.

## 4. Data Source

The data were obtained from the TCGA Breast Invasive Carcinoma study (`brca_tcga`) through cBioPortal.

The analysis used primary tumour samples and selected mutation data for four genes:

| Gene | Entrez Gene ID |
|---|---:|
| TP53 | 7157 |
| BRCA1 | 672 |
| BRCA2 | 675 |
| PIK3CA | 5290 |

The original clinical dataset contained 1,108 records. After restricting the analysis to primary tumour samples, 1,101 samples were retained.

## 5. Tools and Technologies

The project uses:

- R (version 4.6.0)
- RStudio
- cBioPortal
- cBioPortalData
- dplyr
- tidyr
- readr
- ggplot2
- survival
- survminer
- here
- janitor
- Git and GitHub

## 6. Analytical Approach

The analysis was conducted in four main stages.

### 6.1 Data Acquisition

Clinical and mutation data were retrieved from cBioPortal. The selected gene mutation records and clinical information were saved for subsequent analysis.

### 6.2 Data Quality Control and Processing

Quality checks were performed to assess missing values, duplicate records, sample identifiers, and data consistency.

The dataset was restricted to primary tumour samples. Relevant clinical variables, including age, tumour stage, overall survival time, and survival status, were prepared for statistical analysis.

### 6.3 Mutation and Clinical Association Analysis

Mutation frequencies were calculated for the four selected genes.

Chi-square tests and Fisher's exact tests were used, where appropriate, to investigate associations between mutation status and tumour stage.

### 6.4 Survival Analysis

Kaplan-Meier survival analysis and log-rank tests were used to compare overall survival between mutation-status groups.

Cox proportional hazards regression was also performed to examine the relationship between mutation status and survival while accounting for age and tumour stage.

A stage-stratified Cox model was fitted to address the proportional hazards concern identified in the initial model.

## 7. Results and Key Findings

### 7.1 Mutation Frequencies

The mutation analysis produced the following results:

| Gene | Mutated Samples | Frequency |
|---|---:|---:|
| TP53 | 300 | 27.25% |
| BRCA1 | 13 | 1.18% |
| BRCA2 | 15 | 1.36% |
| PIK3CA | 317 | 28.79% |

PIK3CA and TP53 had the highest mutation frequencies among the four genes examined.

The gene-specific mutation counts are not mutually exclusive because a sample may contain mutations in more than one gene.

### 7.2 Association Between Mutation Status and Tumour Stage

No statistically significant association was observed between the combined mutation status of the selected genes and tumour stage.

The combined analysis produced a p-value of 0.478.

The gene-specific analyses also showed no statistically significant associations with tumour stage.

### 7.3 Kaplan-Meier Survival Analysis

The Kaplan-Meier analysis did not identify a statistically significant difference in overall survival between the combined mutation-status groups.

The log-rank test produced:

- Chi-square statistic: 0.607
- Degrees of freedom: 1
- P-value: 0.436
- Number of samples: 1,083

The gene-specific survival comparisons also showed no statistically significant differences.

### 7.4 Cox Proportional Hazards Regression

A multivariable Cox regression model was used to examine the association between mutation status and overall survival while accounting for age and tumour stage.

The initial model showed that age and tumour stage were associated with overall survival. However, the proportional hazards assumption was not satisfied for tumour stage and the overall model.

A stage-stratified Cox model was therefore fitted.

Key findings from the stratified model:

| Variable | Hazard Ratio | 95% Confidence Interval | P-value |
|---|---:|---:|---:|
| Mutation status | 0.81 | 0.57–1.13 | 0.210 |
| Age (per year) | 1.035 | 1.021–1.048 | <0.001 |

The global proportional hazards test for the stratified model produced a p-value of 0.238.

The combined mutation status of the four selected genes was not significantly associated with overall survival after accounting for age and stratifying by tumour stage.

Age showed a statistically significant association with overall survival.

## 8. Analysis Cohorts

Different sample sizes were used across the analyses because of missing clinical information and the requirements of each statistical method.

| Analysis | Sample Size |
|---|---:|
| Primary tumour cohort | 1,101 |
| Mutation-stage association | 1,073 |
| Kaplan-Meier analysis | 1,083 |
| Cox regression analysis | 1,058 |

The Cox regression analysis included 140 deaths.

## 9. Interpretation

The findings describe the mutation patterns of four selected genes in the TCGA breast cancer dataset.

TP53 and PIK3CA had higher mutation frequencies than BRCA1 and BRCA2 in the analysed samples.

However, the statistical analyses did not demonstrate significant associations between the combined mutation status of the selected genes and tumour stage or overall survival.

Age showed a significant association with overall survival in the adjusted Cox model.

These findings should be interpreted within the scope of the selected genes and the available clinical data. The absence of detected mutations in these four genes does not establish that a sample is genomically wild type.

## 10. Limitations

The project has several limitations:

- The analysis focused on only four selected genes.
- Missing clinical information reduced the sample size in some analyses.
- Mutation status was based on the selected genes rather than the entire tumour genome.
- The number of samples with BRCA1 and BRCA2 mutations was relatively small.
- The analysis used observational data and cannot establish causal relationships.
- Findings from the TCGA dataset may not represent all breast cancer populations.

## 11. Project Structure

The repository is organised as follows:

```text
TCGA-BRCA-Bioinformatics/
|
|-- data/
|   |-- raw/
|   |-- processed/
|
|-- scripts/
|   |-- 01_project_setup.R
|   |-- 02_data_acquisition.R
|   |-- 03_data_processing.R
|   |-- 04_statistical_analysis.R
|
|-- results/
|   |-- figures/
|   |-- tables/
|
|-- report/
|
|-- .gitignore
|-- README.md
|-- TCGA-BRCA-Bioinformatics.Rproj
```

## 12. Reproducibility

To reproduce the analysis:

1. Install R and RStudio.
2. Clone or download this repository.
3. Open the RStudio project file.
4. Run the scripts in numerical order, beginning with `01_project_setup.R`.
5. Follow the instructions and package requirements documented in the scripts.

The project uses cBioPortal data. Re-running the acquisition script may require an active internet connection.

Raw and processed datasets are excluded from version control through `.gitignore`.

## 13. Project Status

The data acquisition, processing, quality control, and statistical analysis stages have been completed.

The project is currently at the results interpretation and reporting stage.

## 14. Author

**Clifford Yakubu Biney**

Bioinformatics | Molecular Biology | Cancer Research

GitHub: [Cliffordhub](https://github.com/Cliffordhub)
