
# ============================================================
# SCRIPT 03: DATA PROCESSING
# Project: TCGA-BRCA Bioinformatics
# Purpose: Clean clinical data and prepare mutation data
# ============================================================

# 1. Load required packages ----------------------------------

library(dplyr)
library(here)
library(readr)

# 2. Load the original datasets -------------------------------

clinical_data <- readRDS(
  here("data", "raw", "primary_clinical_data.rds")
)

mutation_data <- readRDS(
  here("data", "raw", "selected_gene_mutations.rds")
)

cat("Original datasets loaded successfully.\n")

# 3. Check required columns -----------------------------------

required_clinical <- c(
  "patientId",
  "sampleId",
  "AGE",
  "OS_MONTHS",
  "OS_STATUS",
  "AJCC_PATHOLOGIC_TUMOR_STAGE"
)

required_mutation <- c(
  "sampleId",
  "entrezGeneId"
)

if (!all(required_clinical %in% names(clinical_data))) {
  stop("Required clinical columns are missing.")
}

if (!all(required_mutation %in% names(mutation_data))) {
  stop("Required mutation columns are missing.")
}

# 4. Process clinical data -----------------------------------

clinical_processed <- clinical_data %>%
  mutate(
    AGE = as.numeric(AGE),
    OS_MONTHS = as.numeric(OS_MONTHS),
    
    # Preserve the original stage information
    stage_original = as.character(
      AJCC_PATHOLOGIC_TUMOR_STAGE
    ),
    
    # Standardise stages into four main groups.
    # Stage X and unrecognised values become NA.
    stage_clean = toupper(trimws(stage_original)),
    
    stage_group = case_when(
      grepl("^STAGE I(A|B)?$", stage_clean) ~ "STAGE I",
      grepl("^STAGE II(A|B)?$", stage_clean) ~ "STAGE II",
      grepl("^STAGE III(A|B|C)?$", stage_clean) ~ "STAGE III",
      grepl("^STAGE IV$", stage_clean) ~ "STAGE IV",
      TRUE ~ NA_character_
    ),
    
    # Convert survival status into an event indicator:
    # 1 = deceased; 0 = living
    survival_status_original = as.character(OS_STATUS),
    
    survival_event = case_when(
      grepl("DECEASED", OS_STATUS, ignore.case = TRUE) ~ 1,
      grepl("LIVING", OS_STATUS, ignore.case = TRUE) ~ 0,
      TRUE ~ NA_real_
    )
  )

# 5. Process mutation data -----------------------------------

gene_mapping <- c(
  "7157" = "TP53",
  "672" = "BRCA1",
  "675" = "BRCA2",
  "5290" = "PIK3CA"
)

mutation_annotated <- mutation_data %>%
  mutate(
    gene_name = unname(
      gene_mapping[as.character(entrezGeneId)]
    )
  ) %>%
  filter(!is.na(gene_name))

# 6. Create sample-level mutation indicators -----------------

mutation_status_by_sample <- mutation_annotated %>%
  distinct(sampleId, gene_name) %>%
  mutate(mutated = 1L) %>%
  tidyr::pivot_wider(
    names_from = gene_name,
    values_from = mutated,
    values_fill = 0
  )

# Ensure all four gene columns exist
for (gene in c("TP53", "BRCA1", "BRCA2", "PIK3CA")) {
  if (!gene %in% names(mutation_status_by_sample)) {
    mutation_status_by_sample[[gene]] <- 0L
  }
}

# 7. Combine clinical and mutation data ----------------------

analysis_data <- clinical_processed %>%
  left_join(
    mutation_status_by_sample,
    by = "sampleId"
  ) %>%
  mutate(
    across(
      all_of(c("TP53", "BRCA1", "BRCA2", "PIK3CA")),
      ~ coalesce(as.integer(.x), 0L)
    ),
    
    any_selected_gene_mutation =
      TP53 == 1 | BRCA1 == 1 | BRCA2 == 1 | PIK3CA == 1,
    
    mutation_status = if_else(
      any_selected_gene_mutation,
      "Mutated",
      "No detected mutation"
    ),
    
    mutation_status = factor(
      mutation_status,
      levels = c("Mutated", "No detected mutation")
    )
  )

# 8. Save processed datasets ---------------------------------

saveRDS(
  clinical_processed,
  here("data", "processed", "clinical_processed.rds")
)

saveRDS(
  mutation_annotated,
  here("data", "processed", "mutation_annotated.rds")
)

saveRDS(
  mutation_status_by_sample,
  here("data", "processed", "mutation_status_by_sample.rds")
)

saveRDS(
  analysis_data,
  here("data", "processed", "analysis_dataset.rds")
)

# 9. Generate mutation summary -------------------------------

mutation_summary <- analysis_data %>%
  summarise(
    total_samples = n(),
    TP53 = sum(TP53),
    BRCA1 = sum(BRCA1),
    BRCA2 = sum(BRCA2),
    PIK3CA = sum(PIK3CA),
    samples_with_selected_gene_mutation =
      sum(any_selected_gene_mutation)
  )

write_csv(
  mutation_summary,
  here("results", "tables", "mutation_status_summary.csv")
)

# 10. Data quality checks ------------------------------------

cat("\n====================================\n")
cat("DATA PROCESSING CHECKS\n")
cat("====================================\n")

cat("Clinical samples:", nrow(clinical_processed), "\n")
cat("Mutation records:", nrow(mutation_annotated), "\n")
cat("Analysis samples:", nrow(analysis_data), "\n")

cat(
  "Duplicate sample IDs:",
  sum(duplicated(analysis_data$sampleId)),
  "\n"
)

cat("\nStage distribution:\n")
print(table(analysis_data$stage_group, useNA = "ifany"))

cat("\nMissing data:\n")

print(
  analysis_data %>%
    summarise(
      missing_age = sum(is.na(AGE)),
      missing_stage = sum(is.na(stage_group)),
      missing_survival_time = sum(is.na(OS_MONTHS)),
      nonpositive_survival_time = sum(
        !is.na(OS_MONTHS) & OS_MONTHS <= 0
      ),
      missing_survival_status = sum(is.na(survival_event))
    ),
  width = Inf
)

cat("\nMutation summary:\n")
print(mutation_summary)

cat("\nDATA PROCESSING COMPLETED\n")
cat(
  "Processed analysis dataset:",
  nrow(analysis_data),
  "samples and",
  ncol(analysis_data),
  "variables.\n"
)
