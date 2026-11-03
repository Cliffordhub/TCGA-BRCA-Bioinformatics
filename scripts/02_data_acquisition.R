
# ============================================================
# 02. DATA ACQUISITION AND QUALITY CONTROL
# TCGA-BRCA Bioinformatics Project
#
# Purpose:
# 1. Load existing datasets where available.
# 2. Retrieve missing datasets from cBioPortal.
# 3. Validate dataset structure and integrity.
# 4. Record data provenance and software versions.
# ============================================================


# 1. Load required packages ----------------------------------

required_packages <- c(
  "cBioPortalData",
  "dplyr",
  "here",
  "readr"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    ". Run the project setup script first."
  )
}

library(cBioPortalData)
library(dplyr)


# 2. Define project settings ---------------------------------

study_id <- "brca_tcga"
mutation_profile_id <- "brca_tcga_mutations"

gene_ids <- c(7157, 672, 675, 5290)

gene_names <- c(
  "TP53",
  "BRCA1",
  "BRCA2",
  "PIK3CA"
)

raw_dir <- here::here("data", "raw")
processed_dir <- here::here("data", "processed")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)


# 3. Define dataset locations ---------------------------------

clinical_file <- file.path(raw_dir, "clinical_data.rds")
primary_file <- file.path(raw_dir, "primary_clinical_data.rds")
mutation_file <- file.path(raw_dir, "selected_gene_mutations.rds")

clinical_exists <- file.exists(clinical_file)
primary_exists <- file.exists(primary_file)
mutation_exists <- file.exists(mutation_file)

message("\nChecking existing datasets...")


# 4. Load existing datasets ----------------------------------

if (clinical_exists) {
  clinical_data <- readRDS(clinical_file)
  message("Loaded existing clinical dataset.")
}

if (primary_exists) {
  primary_data <- readRDS(primary_file)
  message("Loaded existing primary tumour dataset.")
}

if (mutation_exists) {
  mutation_table <- readRDS(mutation_file)
  message("Loaded existing mutation dataset.")
}


# 5. Retrieve missing datasets -------------------------------

needs_clinical_download <- !clinical_exists || !primary_exists
needs_mutation_download <- !mutation_exists

if (needs_clinical_download || needs_mutation_download) {
  
  message("\nSome datasets are missing. Connecting to cBioPortal...")
  
  cbio <- tryCatch(
    cBioPortal(),
    error = function(e) {
      stop("cBioPortal connection failed: ", conditionMessage(e))
    }
  )
  
  if (is.null(cbio)) {
    stop("Could not establish a cBioPortal connection.")
  }
  
  # Retrieve clinical data only when necessary
  if (needs_clinical_download) {
    
    message("Retrieving clinical data...")
    
    clinical_data <- clinicalData(
      cbio,
      studyId = study_id
    )
    
    if (nrow(clinical_data) == 0) {
      stop("Clinical data retrieval returned no records.")
    }
    
    if (!"SAMPLE_TYPE" %in% names(clinical_data)) {
      stop("SAMPLE_TYPE is missing from the clinical dataset.")
    }
    
    primary_data <- clinical_data %>%
      filter(SAMPLE_TYPE == "Primary")
    
    if (nrow(primary_data) == 0) {
      stop("No primary tumour samples were identified.")
    }
    
    saveRDS(clinical_data, clinical_file)
    saveRDS(primary_data, primary_file)
    
    message("Clinical and primary tumour datasets saved.")
  }
  
  # Retrieve mutations only when necessary
  if (needs_mutation_download) {
    
    message("Retrieving selected-gene mutation data...")
    
    if (!exists("primary_data")) {
      primary_data <- readRDS(primary_file)
    }
    
    profiles <- molecularProfiles(
      cbio,
      studyId = study_id
    )
    
    if (!mutation_profile_id %in% profiles$molecularProfileId) {
      stop("The required mutation profile was not found.")
    }
    
    mutation_data <- mutationData(
      api = cbio,
      molecularProfileIds = mutation_profile_id,
      entrezGeneIds = gene_ids,
      sampleIds = unique(primary_data$sampleId)
    )
    
    if (!mutation_profile_id %in% names(mutation_data)) {
      stop("The expected mutation profile was not returned.")
    }
    
    mutation_table <- mutation_data[[mutation_profile_id]]
    
    if (is.null(mutation_table) || nrow(mutation_table) == 0) {
      stop("No mutation records were returned.")
    }
    
    saveRDS(mutation_table, mutation_file)
    
    message("Selected-gene mutation dataset saved.")
  }
  
} else {
  
  message("All required datasets already exist.")
  message("No download is necessary.")
}


# 6. Validate dataset structure ------------------------------

message("\nDATA STRUCTURE CHECKS")

if (!is.data.frame(clinical_data) ||
    !is.data.frame(primary_data) ||
    !is.data.frame(mutation_table)) {
  stop("One or more datasets are not data frames.")
}

required_clinical_columns <- c(
  "patientId",
  "sampleId",
  "SAMPLE_TYPE"
)

required_mutation_columns <- c(
  "sampleId",
  "entrezGeneId"
)

if (!all(required_clinical_columns %in% names(clinical_data))) {
  stop("Required columns are missing from clinical_data.")
}

if (!all(required_clinical_columns %in% names(primary_data))) {
  stop("Required columns are missing from primary_data.")
}

if (!all(required_mutation_columns %in% names(mutation_table))) {
  stop("Required columns are missing from mutation_table.")
}

if (nrow(clinical_data) == 0 ||
    nrow(primary_data) == 0 ||
    nrow(mutation_table) == 0) {
  stop("One or more datasets are empty.")
}

message("Dataset structures are valid.")


# 7. Validate primary tumour samples --------------------------

message("\nPRIMARY TUMOUR QUALITY CHECKS")

if (any(primary_data$SAMPLE_TYPE != "Primary", na.rm = TRUE)) {
  stop("Non-primary samples were found in primary_data.")
}

duplicate_patient_ids <- sum(
  duplicated(primary_data$patientId)
)

duplicate_sample_ids <- sum(
  duplicated(primary_data$sampleId)
)

if (duplicate_patient_ids > 0) {
  stop(
    "Duplicate patient IDs detected: ",
    duplicate_patient_ids
  )
}

if (duplicate_sample_ids > 0) {
  stop(
    "Duplicate sample IDs detected: ",
    duplicate_sample_ids
  )
}

message("All samples are primary tumour samples.")
message("Duplicate patient IDs: ", duplicate_patient_ids)
message("Duplicate sample IDs: ", duplicate_sample_ids)


# 8. Validate mutation records -------------------------------

message("\nMUTATION DATA QUALITY CHECKS")

mutation_table$entrezGeneId <- as.integer(
  as.character(mutation_table$entrezGeneId)
)

if (anyNA(mutation_table$entrezGeneId)) {
  stop("Missing or invalid Entrez Gene IDs in mutation data.")
}

unexpected_genes <- setdiff(
  unique(mutation_table$entrezGeneId),
  gene_ids
)

if (length(unexpected_genes) > 0) {
  stop(
    "Unexpected gene IDs found: ",
    paste(unexpected_genes, collapse = ", ")
  )
}

primary_sample_ids <- unique(primary_data$sampleId)

mutation_sample_ids <- unique(mutation_table$sampleId)

unexpected_samples <- setdiff(
  mutation_sample_ids,
  primary_sample_ids
)

if (length(unexpected_samples) > 0) {
  stop(
    "Mutation records contain samples outside the primary tumour cohort."
  )
}

message("Mutation records contain only the selected genes.")
message("All mutation sample IDs belong to the primary tumour cohort.")


# 9. Record gene-level mutation counts -----------------------

gene_summary <- data.frame(
  gene_name = gene_names,
  entrez_gene_id = gene_ids,
  mutated_samples = vapply(
    gene_ids,
    function(id) {
      length(unique(
        mutation_table$sampleId[
          mutation_table$entrezGeneId == id
        ]
      ))
    },
    integer(1)
  )
)

gene_summary$frequency_percent <- round(
  gene_summary$mutated_samples /
    nrow(primary_data) * 100,
  2
)

message("\nSELECTED GENE MUTATION SUMMARY")
print(gene_summary, row.names = FALSE)


# 10. Record acquisition information -------------------------

acquisition_information <- data.frame(
  study_id = study_id,
  mutation_profile = mutation_profile_id,
  retrieval_or_load_date = as.character(Sys.Date()),
  clinical_records = nrow(clinical_data),
  primary_samples = nrow(primary_data),
  mutation_records = nrow(mutation_table),
  duplicate_patient_ids = duplicate_patient_ids,
  duplicate_sample_ids = duplicate_sample_ids,
  selected_genes = paste(gene_names, collapse = ", ")
)

readr::write_csv(
  acquisition_information,
  file.path(processed_dir, "data_acquisition_summary.csv")
)

readr::write_csv(
  gene_summary,
  file.path(processed_dir, "acquisition_gene_summary.csv")
)


# 11. Record software versions and session information -------

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package) {
      as.character(packageVersion(package))
    },
    character(1)
  )
)

readr::write_csv(
  package_versions,
  file.path(processed_dir, "acquisition_package_versions.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(processed_dir, "acquisition_session_info.txt")
)


# 12. Record data provenance ---------------------------------

provenance <- c(
  "PROJECT: TCGA-BRCA Bioinformatics",
  "DATA SOURCE: cBioPortal",
  paste("STUDY ID:", study_id),
  paste("MUTATION PROFILE:", mutation_profile_id),
  paste("SELECTED GENES:", paste(gene_names, collapse = ", ")),
  paste("ACQUISITION/LOAD DATE:", Sys.Date()),
  paste("CLINICAL RECORDS:", nrow(clinical_data)),
  paste("PRIMARY TUMOUR SAMPLES:", nrow(primary_data)),
  paste("MUTATION RECORDS:", nrow(mutation_table)),
  "DATA TYPE: Public research data"
)

writeLines(
  provenance,
  file.path(processed_dir, "data_provenance.txt")
)


# 13. Final quality checks -----------------------------------

stopifnot(
  nrow(clinical_data) == 1108,
  nrow(primary_data) == 1101,
  nrow(mutation_table) == 686,
  duplicate_patient_ids == 0,
  duplicate_sample_ids == 0,
  all(mutation_table$entrezGeneId %in% gene_ids),
  all(mutation_table$sampleId %in% primary_sample_ids)
)

message("\n==========================================")
message("DATA ACQUISITION AND QC COMPLETED")
message("==========================================")

message("Clinical records: ", nrow(clinical_data))
message("Primary tumour samples: ", nrow(primary_data))
message("Mutation records: ", nrow(mutation_table))
message("Duplicate patient IDs: ", duplicate_patient_ids)
message("Duplicate sample IDs: ", duplicate_sample_ids)

message("\nSelected genes:")
message(paste(gene_names, collapse = ", "))

message("\nAll acquisition quality checks passed.")
message("==========================================")
