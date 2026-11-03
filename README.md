#### **Exploring Molecular Alterations and Their Association with Clinical Characteristics and Overall Survival in Breast Cancer Using TCGA Data**

#### 

###### **Project Overview

###### This project explores selected genetic alterations in breast cancer using publicly available data from The Cancer Genome Atlas (TCGA). It focuses on four genes: TP53, BRCA1, BRCA2, and PIK3CA, examining their mutation frequencies, distribution across tumour stages, and potential association with overall survival.

###### 

###### The analysis was conducted using R and data accessed through cBioPortal.



##### **Research Question**

###### Are mutations in selected breast cancer-associated genes related to clinical characteristics and overall survival among patients with breast cancer?



##### **Study Objectives**

* ###### Determine the mutation frequencies of TP53, BRCA1, BRCA2, and PIK3CA.
* ###### Describe the distribution of mutation types among the selected genes.
* ###### Explore the relationship between mutation status and tumour stage.
* ###### Examine overall survival patterns according to mutation status.





##### **Data Source and Retrieval**



###### This project uses publicly available data from the TCGA Breast Invasive Carcinoma (TCGA-BRCA) study, accessed through cBioPortal(https://www.cbioportal.org/).



###### \- Dataset: TCGA Breast Invasive Carcinoma

###### \- Study ID: brca\_tcga

###### \- Data platform: cBioPortal

###### \- Samples: 1,101 primary tumour samples used for the clinical and mutation analysis.





##### **Reproducing the Data Retrieval**

###### The data retrieval and preparation procedures are documented in the R scripts provided in the scripts/ directory.



###### To reproduce the analysis:



###### 1\. Install R and RStudio.

###### 2\. Clone or download this repository.

###### 3\. Open the TCGA-BRCA-Bioinformatics.Rproj file in RStudio.

###### 4\. Run scripts/01\_project\_setup.R to prepare the working environment and install the required packages.

###### 5\. Run scripts/02\_data\_acquisition.R to retrieve and prepare the study data through cBioPortal.



###### For direct access to the original dataset, visit cBioPortal (https://www.cbioportal.org/) and search for **TCGA Breast Invasive Carcinoma.**

###### 

###### Note: The original study data are not stored in this repository. They are retrieved through the documented workflow. The number of samples available for individual analyses may vary because of missing data and filtering criteria.





##### &#x20; **Selected Genes**

###### &#x20;  **Gene**	       **Biological role**



1. ###### TP53	        Tumour suppressor involved in DNA damage response
2. ###### BRCA1	DNA repair through homologous recombination
3. ###### BRCA2	DNA repair through homologous recombination
4. ###### PIK3CA	Regulates PI3K/AKT cell growth signalling





##### **Tools and Technologies**

###### \- R and RStudio

###### \- cBioPortal

###### \- tidyverse

###### \- survival

###### \- survminer

###### \- janitor

###### \- ggplot2





##### **Analytical Approach**



###### 1\. Data acquisition from cBioPortal.

###### 2\. Data cleaning and preparation.

###### 3\. Descriptive analysis of mutation frequencies and mutation types.

###### 4\. Chi-square and other appropriate statistical tests to explore associations with tumour stage.

###### 5\. Kaplan–Meier survival analysis and log-rank testing.





##### **Preliminary Findings**



* ###### PIK3CA and TP53 had the highest observed mutation frequencies among the four selected genes.
* ###### BRCA1 and BRCA2 mutations were less frequent in the analysed samples.
* ###### No statistically significant association was observed between the combined mutation status of the four genes and tumour stage (p = 0.4775).
* ###### Kaplan–Meier analysis did not identify a statistically significant difference in overall survival between the combined mutation-status groups (log-rank p = 0.44).





###### Combined mutation status was defined using the presence or absence of detected mutations in the four selected genes. The absence of detected mutations in these genes does not imply that a sample was free of other genomic alterations.

###### 

##### **Project Structure**



###### TCGA-BRCA-Bioinformatics/

###### ├── data/

###### │   ├── raw/

###### │   └── processed/

###### ├── scripts/

###### │   ├── 01\_project\_setup.R

###### │   └── 02\_data\_acquisition.R

###### ├── results/

###### │   ├── figures/

###### │   └── tables/

###### ├── report/

###### ├── README.md

###### ├── .gitignore

###### └── TCGA-BRCA-Bioinformatics.Rproj



###### **Limitations**



* ###### The analysis focuses on four selected genes and does not represent all genomic alterations in breast cancer.
* ###### BRCA1 and BRCA2 had relatively small mutation groups, limiting the strength of gene-specific comparisons.
* ###### Missing clinical information may affect some analyses.
* ###### The analysis is observational and cannot establish causation.
* ###### Mutation detection does not automatically establish that a variant is pathogenic or functionally damaging.





###### **Project Status**

###### 

###### This project is part of an ongoing learning exercise in cancer bioinformatics, R programming, statistical analysis, and the interpretation of public cancer genomics datasets.

###### 

###### **Author**

###### Clifford Yakubu Biney

###### 

###### Bioinformatics | Molecular Biology | Cancer Genomics

