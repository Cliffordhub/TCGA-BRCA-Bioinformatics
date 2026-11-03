
# ==========================================
# 02. CONNECTING TO CBIOPORTAL
# TCGA-BRCA Bioinformatics Project
# ==========================================

# Install Bioconductor manager if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install cBioPortalData
if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
  BiocManager::install("cBioPortalData")
}

# Load the package
library(cBioPortalData)

# Connect to cBioPortal
cbio <- cBioPortal()

# Retrieve available cancer studies
studies <- getStudies(cbio)

# Display the first six studies
head(studies)

# Search for breast cancer studies

library(dplyr)

brca_studies <- studies %>%
  filter(
    grepl(
      "breast|brca",
      paste(name, description, studyId),
      ignore.case = TRUE
    )
  )

# Display relevant studies
print(brca_studies)
print(n= ...)
write.csv(
  brca_studies,
  here::here(
    "data",
    "processed",
    "available_breast_cancer_studies.csv"
  ),
  row.names = FALSE
)

# ==========================================
# 03. EXAMINING BREAST CANCER STUDIES
# ==========================================

# Display the structure of the study dataset
str(brca_studies)

# Display column names
colnames(brca_studies)

# Display the first 10 records
head(brca_studies, 10)

# Count the number of available studies
nrow(brca_studies)



# Search specifically for TCGA breast cancer studies

tcga_brca <- brca_studies %>%
  filter(
    grepl("tcga", studyId, ignore.case = TRUE)
  )

# Display important study information
tcga_brca %>%
  select(any_of(c(
    "studyId",
    "name",
    "description",
    "allSampleCount",
    "publicStudy",
    "readPermission"
  ))) %>%
  print(n = Inf)

# Set the selected study ID
study_id <- "brca_tcga"

# Examine available molecular profiles
profiles <- molecularProfiles(
  cbio,
  studyId = study_id
)

# Display the available profiles
print(profiles)

# Examine available clinical information
clinical_data <- clinicalData(
  cbio,
  studyId = study_id
)

# Examine the clinical dataset
dim(clinical_data)
head(clinical_data)



# ==========================================
# 04. CLINICAL DATA EXPLORATION
# ==========================================

# 1. Check the number of rows and columns
dim(clinical_data)

# 2. Display all column names
colnames(clinical_data)

# 3. Check the first few records
head(clinical_data)

# 4. Check for missing values
missing_values <- colSums(is.na(clinical_data))

# Display variables with missing values
sort(missing_values, decreasing = TRUE)[1:15]

# 5. Check duplicate patient IDs
sum(duplicated(clinical_data$patientId))

# 6. Count unique patient IDs
length(unique(clinical_data$patientId))

# 7. Examine cancer stages
table(
  clinical_data$AJCC_PATHOLOGIC_TUMOR_STAGE,
  useNA = "ifany"
)

# 8. Examine patient ages
summary(as.numeric(clinical_data$AGE))




# 9. Identify duplicated patient records

duplicated_patients <- clinical_data %>%
  filter(duplicated(patientId) | duplicated(patientId, fromLast = TRUE)) %>%
  arrange(patientId)

# Display duplicated records
print(duplicated_patients)

# 10. Examine patient and sample identifiers

duplicated_patients %>%
  select(patientId, sampleId, SAMPLE_TYPE, AGE,
         AJCC_PATHOLOGIC_TUMOR_STAGE)

# 11. Examine overall survival information

# Check the first few survival records
head(clinical_data %>%
       select(patientId, OS_MONTHS, OS_STATUS))

# Count missing overall survival values
colSums(is.na(clinical_data[, c("OS_MONTHS", "OS_STATUS")]))

# Examine the different survival outcome categories
table(clinical_data$OS_STATUS, useNA = "ifany")

# Summarise overall survival time
summary(as.numeric(clinical_data$OS_MONTHS))



# 12. Investigate records with negative overall survival time

clinical_data %>%
  mutate(OS_MONTHS_NUM = as.numeric(OS_MONTHS)) %>%
  filter(OS_MONTHS_NUM < 0) %>%
  select(patientId, sampleId, SAMPLE_TYPE, OS_MONTHS, OS_STATUS)

# 13. Examine survival information for patients
# with both primary and metastatic samples

clinical_data %>%
  filter(patientId %in% duplicated_patients$patientId) %>%
  select(patientId, sampleId, SAMPLE_TYPE,
         OS_MONTHS, OS_STATUS) %>%
  arrange(patientId, sampleId)

# 14. Examine the distribution of sample types

table(clinical_data$SAMPLE_TYPE, useNA = "ifany")

# 15. Check whether the primary tumour records
# have unique patient IDs

primary_data <- clinical_data %>%
  filter(SAMPLE_TYPE == "Primary")

dim(primary_data)

sum(duplicated(primary_data$patientId))

length(unique(primary_data$patientId))



# 16. Create a patient-level clinical dataset
# using primary tumour samples

clinical_primary <- primary_data

# 17. Convert age and survival time to numeric values

clinical_primary <- clinical_primary %>%
  mutate(
    AGE = as.numeric(AGE),
    OS_MONTHS = as.numeric(OS_MONTHS)
  )

# 18. Check the new dataset

dim(clinical_primary)

# Confirm one row per patient
sum(duplicated(clinical_primary$patientId))

# Check missing values in important variables
colSums(is.na(
  clinical_primary[, c(
    "AGE",
    "AJCC_PATHOLOGIC_TUMOR_STAGE",
    "OS_MONTHS",
    "OS_STATUS"
  )]
))

# 19. Identify records with invalid survival time

clinical_primary %>%
  filter(OS_MONTHS < 0) %>%
  select(patientId, OS_MONTHS, OS_STATUS)



# 20. Inspect the molecular profiles available
# in the selected TCGA-BRCA study

print(profiles)

# 21. Inspect the study ID being used

study_id

# 22. Examine the available profile identifiers

profiles$ molecularProfileId



# 23. Select the mutation profile

mutation_profile <- "brca_tcga_mutations"

# 24. Retrieve mutation data from cBioPortal

mutation_data <- molecularData(
  cbio,
  studyId = study_id,
  molecularProfileId = mutation_profile
)

# 25. Inspect the retrieved data

dim(mutation_data)

head(mutation_data)

colnames(mutation_data)


# Check the arguments accepted by molecularData()
args(molecularData)

# Check the installed package version
packageVersion("cBioPortalData")

# Check the function documentation
?molecularData



# 23. Select the mutation profile

mutation_profile <- "brca_tcga_mutations"

# 24. Retrieve mutation data

mutation_data <- molecularData(
  api = cbio,
  molecularProfileIds = mutation_profile
)

# 25. Inspect the retrieved data

dim(mutation_data)

head(mutation_data)

colnames(mutation_data)




# 23. Select the mutation profile
mutation_profile <- "brca_tcga_mutations"

# 24. Define Entrez Gene IDs for our four genes
gene_ids <- c("7157", "672", "675", "5290")

# 25. Retrieve mutation data for the selected genes
mutation_data <- molecularData(
  api = cbio,
  molecularProfileIds = mutation_profile,
  entrezGeneIds = gene_ids
)

# 26. Inspect the retrieved data
dim(mutation_data)
head(mutation_data)
colnames(mutation_data)



# 23. Select the mutation profile
mutation_profile <- "brca_tcga_mutations"

# 24. Define Entrez Gene IDs
gene_ids <- c("7157", "672", "675", "5290")

# 25. Extract primary tumour sample IDs
primary_sample_ids <- unique(primary_data$sampleId)

# Check the number of samples
length(primary_sample_ids)

# 26. Retrieve mutation data
mutation_data <- molecularData(
  api = cbio,
  molecularProfileIds = mutation_profile,
  entrezGeneIds = gene_ids,
  sampleIds = primary_sample_ids
)

# 27. Inspect the retrieved data
dim(mutation_data)
head(mutation_data)
colnames(mutation_data)


# 1. Confirm the study ID
study_id

# 2. Display the complete molecular profile information
print(profiles)

# 3. Retrieve the profiles again for the selected study
profiles_check <- molecularProfiles(
  cbio,
  studyId = study_id
)

print(profiles_check)

# 4. Check the profile IDs returned
profiles_check$molecularProfileId



# 28. Retrieve mutations for a small test group

test_sample_ids <- primary_sample_ids[1:5]

test_mutations <- mutationData(
  api = cbio,
  molecularProfileIds = "brca_tcga_mutations",
  entrezGeneIds = c(7157, 672, 675, 5290),
  sampleIds = test_sample_ids
)

# 29. Inspect the results

dim(test_mutations)

head(test_mutations)

colnames(test_mutations)


# Extract the tibble from the list
test_mutation_table <- test_mutations$brca_tcga_mutations

# Inspect its dimensions
dim(test_mutation_table)

# Inspect its column names
colnames(test_mutation_table)

# View all three mutation records
View(test_mutation_table)

test_mutation_table %>%
  count(entrezGeneId)

test_mutation_table %>%
  count(entrezGeneId)

class(test_mutations)

names(test_mutations)

class(test_mutations$brca_tcga_mutations)


test_mutation_table <- test_mutations[["brca_tcga_mutations"]]

class(test_mutation_table)

dim(test_mutation_table)

dplyr::count(test_mutation_table, entrezGeneId)


# Retrieve mutations across all primary tumour samples
mutation_data <- mutationData(
  api = cbio,
  molecularProfileIds = "brca_tcga_mutations",
  entrezGeneIds = c(7157, 672, 675, 5290),
  sampleIds = primary_sample_ids
)

# Extract the mutation table
mutation_table <- mutation_data[["brca_tcga_mutations"]]

# Check the results
dim(mutation_table)
head(mutation_table)

library(dplyr)
library(ggplot2)

# Count unique samples with at least one mutation in each gene
gene_mutations <- mutation_table %>%
  group_by(entrezGeneId) %>%
  summarise(
    mutated_samples = n_distinct(sampleId),
    mutation_records = n(),
    .groups = "drop"
  )

# Create a table containing all four genes
gene_summary <- tibble(
  entrezGeneId = c(7157, 672, 675, 5290),
  gene = c("TP53", "BRCA1", "BRCA2", "PIK3CA")
) %>%
  left_join(gene_mutations, by = "entrezGeneId") %>%
  mutate(
    mutated_samples = coalesce(mutated_samples, 0L),
    mutation_records = coalesce(mutation_records, 0L),
    total_samples = length(primary_sample_ids),
    mutation_frequency = round(
      mutated_samples / total_samples * 100, 2
    )
  )

# Display the results
gene_summary

ggplot(gene_summary,
       aes(x = reorder(gene, mutation_frequency),
           y = mutation_frequency)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Mutation Frequency of Selected Genes in TCGA-BRCA",
    x = "Gene",
    y = "Samples with Mutations (%)"
  ) +
  theme_minimal()

ggsave(
  "results/figures/gene_mutation_frequency.png",
  width = 8,
  height = 5,
  dpi = 300
)

ggplot(gene_summary,
       aes(x = reorder(gene, mutation_frequency),
           y = mutation_frequency)) +
  geom_col(fill = "steelblue") +
  geom_text(
    aes(label = paste0(mutation_frequency, "%")),
    hjust = -0.2,
    size = 4
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 35),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Mutation Frequency of Selected Genes in Breast Cancer",
    subtitle = "TCGA-BRCA | 1,101 primary tumour samples",
    x = "Gene",
    y = "Samples with at least one detected mutation (%)"
  ) +
  theme_minimal(base_size = 12)

mutation_table %>%
  count(mutationType, sort = TRUE)

class(mutation_table)

names(mutation_table)

dplyr::count(
  mutation_table,
  mutationType,
  sort = TRUE
)

# Identify unique samples with at least one mutation
mutated_samples <- mutation_table %>%
  distinct(sampleId) %>%
  mutate(mutation_status = "Mutated")

# Check the number of unique mutated samples
nrow(mutated_samples)

clinical_analysis <- clinical_primary %>%
  left_join(mutated_samples, by = "sampleId") %>%
  mutate(
    mutation_status = coalesce(
      mutation_status,
      "No detected mutation"
    )
  )

# Check dimensions
dim(clinical_analysis)

# Count patients in each mutation group
clinical_analysis %>%
  count(mutation_status)

# Check tumour stage distribution
clinical_analysis %>%
  count(AJCC_PATHOLOGIC_TUMOR_STAGE, mutation_status)

class(clinical_analysis)

dim(clinical_analysis)

names(clinical_analysis)

class(clinical_analysis$mutation_status)

head(clinical_analysis$mutation_status)


dplyr::count(clinical_analysis, mutation_status)


dplyr::count(
  clinical_analysis,
  AJCC_PATHOLOGIC_TUMOR_STAGE,
  mutation_status
)


sum(is.na(clinical_analysis$AJCC_PATHOLOGIC_TUMOR_STAGE))


dplyr::count(
  clinical_analysis,
  AJCC_PATHOLOGIC_TUMOR_STAGE,
  mutation_status
) %>%
  print(n = Inf)

clinical_stage <- clinical_analysis %>%
  mutate(
    stage_group = case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage III", "Stage IIIA", "Stage IIIB",
          "Stage IIIC") ~ "Stage III",
      
      AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV" ~ "Stage IV",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(stage_group))


stage_comparison <- clinical_stage %>%
  group_by(stage_group, mutation_status) %>%
  summarise(
    patients = n(),
    .groups = "drop"
  )

stage_comparison


stage_percentages <- clinical_stage %>%
  group_by(stage_group) %>%
  count(mutation_status) %>%
  mutate(
    percentage = round(n / sum(n) * 100, 2)
  )

stage_percentages


stage_percentages <- clinical_stage %>%
  dplyr::group_by(stage_group) %>%
  dplyr::count(mutation_status) %>%
  dplyr::mutate(
    percentage = round(n / sum(n) * 100, 2)
  )

stage_percentages

nrow(clinical_stage)


# Create a contingency table
stage_table <- table(
  clinical_stage$stage_group,
  clinical_stage$mutation_status
)

# View the table
stage_table

# Perform Chi-square test
stage_test <- chisq.test(stage_table)

# View the results
stage_test


# Check whether the Chi-square test assumptions are satisfied
stage_test$expected




# Prepare data for overall survival analysis

survival_data <- clinical_analysis %>%
  mutate(
    OS_MONTHS = as.numeric(OS_MONTHS),
    
    event = case_when(
      OS_STATUS == "1:DECEASED" ~ 1,
      OS_STATUS == "0:LIVING" ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(
    !is.na(OS_MONTHS),
    OS_MONTHS >= 0,
    !is.na(event)
  )

# Check the dimensions of the prepared dataset
dim(survival_data)

# Check the number of patients in each mutation group
dplyr::count(survival_data, mutation_status)

# Check the number of deaths and censored observations
dplyr::count(survival_data, event)

# Summarise the survival time
summary(survival_data$OS_MONTHS)



# Load survival analysis packages
library(survival)
library(survminer)

# Fit Kaplan-Meier survival model
km_fit <- survfit(
  Surv(OS_MONTHS, event) ~ mutation_status,
  data = survival_data
)

# View survival model results
summary(km_fit)



# Plot Kaplan-Meier survival curves

km_plot <- ggsurvplot(
  km_fit,
  data = survival_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  title = "Overall Survival by Mutation Status",
  legend.title = "Mutation Status",
  legend.labs = c("Mutated", "No detected mutation"),
  surv.median.line = "hv"
)

# Display the plot
print(km_plot)

# Save the plot
ggsave(
  "results/figures/kaplan_meier_overall_survival.png",
  print(km_plot),
  width = 9,
  height = 7,
  dpi = 300
)



# Perform the log-rank test
logrank_test <- survdiff(
  Surv(OS_MONTHS, event) ~ mutation_status,
  data = survival_data
)

logrank_test

# Extract median survival estimates
summary(km_fit)$table



# Identify samples with TP53 mutations
tp53_samples <- mutation_table %>%
  filter(entrezGeneId == 7157) %>%
  distinct(sampleId) %>%
  mutate(tp53_status = "Mutated")

# Add TP53 status to the clinical survival data
tp53_survival <- survival_data %>%
  left_join(tp53_samples, by = "sampleId") %>%
  mutate(
    tp53_status = coalesce(
      tp53_status,
      "No detected mutation"
    )
  )

# Check the two groups
dplyr::count(tp53_survival, tp53_status)

# Fit Kaplan-Meier model
tp53_fit <- survfit(
  Surv(OS_MONTHS, event) ~ tp53_status,
  data = tp53_survival
)

# Plot survival curves
ggsurvplot(
  tp53_fit,
  data = tp53_survival,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by TP53 Mutation Status",
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)



# Identify samples with BRCA1 mutations
brca1_samples <- mutation_table %>%
  filter(entrezGeneId == 672) %>%
  distinct(sampleId) %>%
  mutate(brca1_status = "Mutated")

# Add BRCA1 status to survival data
brca1_survival <- survival_data %>%
  left_join(brca1_samples, by = "sampleId") %>%
  mutate(
    brca1_status = coalesce(
      brca1_status,
      "No detected mutation"
    )
  )

# Check group sizes
dplyr::count(brca1_survival, brca1_status)

# Fit Kaplan-Meier model
brca1_fit <- survfit(
  Surv(OS_MONTHS, event) ~ brca1_status,
  data = brca1_survival
)

# Plot survival curves
ggsurvplot(
  brca1_fit,
  data = brca1_survival,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by BRCA1 Mutation Status",
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)



# 1. Identify samples with BRCA1 mutations
brca1_samples <- mutation_table %>%
  filter(entrezGeneId == 672) %>%
  distinct(sampleId) %>%
  mutate(brca1_status = "Mutated")

# 2. Add BRCA1 status to survival data
brca1_survival <- survival_data %>%
  left_join(brca1_samples, by = "sampleId") %>%
  mutate(
    brca1_status = coalesce(
      brca1_status,
      "No detected mutation"
    )
  )

# 3. Check group sizes
dplyr::count(brca1_survival, brca1_status)

# 4. Fit Kaplan-Meier model
brca1_fit <- survfit(
  Surv(OS_MONTHS, event) ~ brca1_status,
  data = brca1_survival
)

# 5. Plot survival curves
ggsurvplot(
  brca1_fit,
  data = brca1_survival,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by BRCA1 Mutation Status",
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)


exists("mutation_table")
exists("survival_data")
exists("brca1_samples")


brca1_samples <- mutation_table %>%
  dplyr::filter(entrezGeneId == 672) %>%
  dplyr::distinct(sampleId) %>%
  dplyr::mutate(brca1_status = "Mutated")

brca1_survival <- survival_data %>%
  dplyr::left_join(brca1_samples, by = "sampleId") %>%
  dplyr::mutate(
    brca1_status = dplyr::coalesce(
      brca1_status,
      "No detected mutation"
    )
  )

exists("brca1_survival")
dplyr::count(brca1_survival, brca1_status)

brca1_fit <- survfit(
  Surv(OS_MONTHS, event) ~ brca1_status,
  data = brca1_survival
)

ggsurvplot(
  brca1_fit,
  data = brca1_survival,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by BRCA1 Mutation Status",
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)



# Identify BRCA2-mutated samples
brca2_samples <- mutation_table %>%
  dplyr::filter(entrezGeneId == 675) %>%
  dplyr::distinct(sampleId) %>%
  dplyr::mutate(brca2_status = "Mutated")

# Add BRCA2 status to survival data
brca2_survival <- survival_data %>%
  dplyr::left_join(brca2_samples, by = "sampleId") %>%
  dplyr::mutate(
    brca2_status = dplyr::coalesce(
      brca2_status,
      "No detected mutation"
    )
  )

# Check group sizes
dplyr::count(brca2_survival, brca2_status)

# Fit Kaplan-Meier model
brca2_fit <- survfit(
  Surv(OS_MONTHS, event) ~ brca2_status,
  data = brca2_survival
)

# Plot survival curves
ggsurvplot(
  brca2_fit,
  data = brca2_survival,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by BRCA2 Mutation Status",
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)



# Identify PIK3CA-mutated samples
pik3ca_samples <- mutation_table %>%
  dplyr::filter(entrezGeneId == 5290) %>%
  dplyr::distinct(sampleId) %>%
  dplyr::mutate(pik3ca_status = "Mutated")

# Add PIK3CA status to survival data
pik3ca_survival <- survival_data %>%
  dplyr::left_join(pik3ca_samples, by = "sampleId") %>%
  dplyr::mutate(
    pik3ca_status = dplyr::coalesce(
      pik3ca_status,
      "No detected mutation"
    )
  )

# Check group sizes
dplyr::count(pik3ca_survival, pik3ca_status)

# Fit Kaplan-Meier model
pik3ca_fit <- survfit(
  Surv(OS_MONTHS, event) ~ pik3ca_status,
  data = pik3ca_survival
)

# Plot survival curves
ggsurvplot(
  pik3ca_fit,
  data = pik3ca_survival,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  title = "Overall Survival by PIK3CA Mutation Status",
  xlab = "Time (months)",
  ylab = "Overall survival probability"
)



# Check available clinical variables
names(survival_data)[
  names(survival_data) %in%
    c("AGE", "AJCC_PATHOLOGIC_TUMOR_STAGE",
      "OS_MONTHS", "event", "mutation_status")
]

# Check missing values
colSums(is.na(
  survival_data[, c(
    "AGE",
    "AJCC_PATHOLOGIC_TUMOR_STAGE",
    "mutation_status"
  )]
))

# Check age distribution
summary(survival_data$AGE)

# Check tumour stage distribution
dplyr::count(
  survival_data,
  AJCC_PATHOLOGIC_TUMOR_STAGE
)




# Prepare tumour stage groups
cox_data <- survival_data %>%
  mutate(
    AGE = as.numeric(AGE),
    
    stage_group = case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage III", "Stage IIIA", "Stage IIIB",
          "Stage IIIC") ~ "Stage III",
      
      AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV" ~ "Stage IV",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(AGE),
    !is.na(stage_group)
  )

# Set Stage I as the reference category
cox_data$stage_group <- factor(
  cox_data$stage_group,
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
)

# Check the final dataset
dim(cox_data)

# Check stage distribution
dplyr::count(cox_data, stage_group)

# Check age summary
summary(cox_data$AGE)


# Fit adjusted Cox proportional hazards model
cox_model <- coxph(
  Surv(OS_MONTHS, event) ~
    mutation_status + AGE + stage_group,
  data = cox_data
)

# Display model results
summary(cox_model)




ph_test <- cox.zph(cox_model)

ph_test
plot(ph_test)


# Plot proportional hazards diagnostics for all variables
plot(ph_test, var = "stage_group")


# Fit stratified Cox regression model
cox_stratified <- coxph(
  Surv(OS_MONTHS, event) ~
    mutation_status + AGE + strata(stage_group),
  data = cox_data
)

# Display results
summary(cox_stratified)

# Check proportional hazards assumption again
ph_stratified <- cox.zph(cox_stratified)

ph_stratified


summary(cox_stratified)


# Identify samples with TP53 mutations
tp53_samples <- mutation_table %>%
  dplyr::filter(entrezGeneId == 7157) %>%
  dplyr::distinct(sampleId) %>%
  dplyr::mutate(tp53_status = "Mutated")

# Combine TP53 status with clinical data
tp53_clinical <- clinical_analysis %>%
  dplyr::left_join(tp53_samples, by = "sampleId") %>%
  dplyr::mutate(
    tp53_status = dplyr::coalesce(
      tp53_status,
      "No detected mutation"
    )
  )

# Group tumour stages
tp53_clinical <- tp53_clinical %>%
  dplyr::mutate(
    stage_group = dplyr::case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage III", "Stage IIIA", "Stage IIIB",
          "Stage IIIC") ~ "Stage III",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV" ~ "Stage IV",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(stage_group))

# Create contingency table
tp53_stage_table <- table(
  tp53_clinical$stage_group,
  tp53_clinical$tp53_status
)

tp53_stage_table

# Chi-square test
chisq.test(tp53_stage_table)



# Identify BRCA1-mutated samples
brca1_samples <- mutation_table %>%
  dplyr::filter(entrezGeneId == 672) %>%
  dplyr::distinct(sampleId) %>%
  dplyr::mutate(brca1_status = "Mutated")

# Combine with clinical data
brca1_clinical <- clinical_analysis %>%
  dplyr::left_join(brca1_samples, by = "sampleId") %>%
  dplyr::mutate(
    brca1_status = dplyr::coalesce(
      brca1_status,
      "No detected mutation"
    ),
    stage_group = dplyr::case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in%
        c("Stage III", "Stage IIIA", "Stage IIIB",
          "Stage IIIC") ~ "Stage III",
      AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV" ~ "Stage IV",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(stage_group))

# Contingency table
brca1_stage_table <- table(
  brca1_clinical$stage_group,
  brca1_clinical$brca1_status
)

brca1_stage_table

# Chi-square test
chisq.test(brca1_stage_table)
