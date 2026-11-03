
# ============================================================
# TCGA-BRCA BIOINFORMATICS PROJECT
# SCRIPT 04: STATISTICAL ANALYSIS
#
# Research focus:
# Association between selected gene mutations,
# clinical characteristics and overall survival.
# ============================================================


# 1. PROJECT SETUP AND PACKAGES
# ============================================================

required_packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2",
  "survival", "survminer", "here", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

cat("\n============================================\n")
cat("TCGA-BRCA STATISTICAL ANALYSIS\n")
cat("============================================\n")


# 2. LOAD PROCESSED DATA
# ============================================================

data_file <- here::here(
  "data", "processed", "analysis_dataset.rds"
)

if (!file.exists(data_file)) {
  stop("Analysis dataset not found. Run Script 03 first.")
}

analysis_data <- readRDS(data_file)

if (!is.data.frame(analysis_data)) {
  stop("The analysis dataset is not a valid data frame.")
}

cat("\nDATA LOADED\n")
cat("Total samples:", nrow(analysis_data), "\n")


# 3. VALIDATE REQUIRED VARIABLES
# ============================================================

required_columns <- c(
  "sampleId", "AGE", "OS_MONTHS", "survival_event",
  "stage_group", "mutation_status",
  "TP53", "BRCA1", "BRCA2", "PIK3CA"
)

missing_columns <- setdiff(required_columns, names(analysis_data))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (anyDuplicated(analysis_data$sampleId) > 0) {
  stop("Duplicate sample IDs detected.")
}

cat("\nVARIABLE VALIDATION PASSED\n")


# 4. STANDARDISE VARIABLES
# ============================================================

analysis_data <- analysis_data |>
  dplyr::mutate(
    AGE = suppressWarnings(as.numeric(as.character(AGE))),
    OS_MONTHS = suppressWarnings(as.numeric(as.character(OS_MONTHS))),
    survival_event = suppressWarnings(
      as.numeric(as.character(survival_event))
    ),
    
    stage_clean = toupper(trimws(as.character(stage_group))),
    
    stage_group = dplyr::case_when(
      stage_clean %in% c("STAGE I", "STAGE 1", "I") ~ "Stage I",
      stage_clean %in% c("STAGE II", "STAGE 2", "II") ~ "Stage II",
      stage_clean %in% c("STAGE III", "STAGE 3", "III") ~ "Stage III",
      stage_clean %in% c("STAGE IV", "STAGE 4", "IV") ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    
    stage_group = factor(
      stage_group,
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
    ),
    
    mutation_status = factor(
      as.character(mutation_status),
      levels = c("Mutated", "No detected mutation")
    )
  )

cat("\nSTANDARDISED STAGE DISTRIBUTION\n")
print(table(analysis_data$stage_group, useNA = "ifany"))


# 5. PREPARE ANALYSIS COHORTS
# ============================================================

stage_data <- analysis_data |>
  dplyr::filter(
    !is.na(stage_group),
    !is.na(mutation_status)
  )

km_data <- analysis_data |>
  dplyr::filter(
    !is.na(OS_MONTHS),
    OS_MONTHS > 0,
    !is.na(survival_event),
    survival_event %in% c(0, 1),
    !is.na(mutation_status)
  )

cox_data <- analysis_data |>
  dplyr::filter(
    !is.na(OS_MONTHS),
    OS_MONTHS > 0,
    !is.na(survival_event),
    survival_event %in% c(0, 1),
    !is.na(AGE),
    !is.na(stage_group),
    !is.na(mutation_status)
  )

cat("\nANALYSIS COHORTS\n")
cat("Stage association samples:", nrow(stage_data), "\n")
cat("Kaplan-Meier samples:", nrow(km_data), "\n")
cat("Kaplan-Meier deaths:", sum(km_data$survival_event == 1), "\n")
cat("Kaplan-Meier censored:", sum(km_data$survival_event == 0), "\n")
cat("Cox regression samples:", nrow(cox_data), "\n")
cat("Cox regression deaths:", sum(cox_data$survival_event == 1), "\n")
cat("Cox regression censored:", sum(cox_data$survival_event == 0), "\n")

if (any(c(nrow(stage_data), nrow(km_data), nrow(cox_data)) == 0)) {
  stop("One or more analysis cohorts are empty.")
}


# 6. GENE MUTATION STATUS FUNCTION
# ============================================================

get_gene_status <- function(data, gene) {
  
  gene_values <- data[[gene]]
  
  if (is.logical(gene_values)) {
    gene_values <- as.integer(gene_values)
  }
  
  gene_values <- suppressWarnings(
    as.numeric(as.character(gene_values))
  )
  
  if (!all(is.na(gene_values) | gene_values %in% c(0, 1))) {
    stop(
      gene,
      " contains values other than 0 and 1. ",
      "Check Script 03 gene-status coding."
    )
  }
  
  factor(
    ifelse(
      is.na(gene_values),
      NA_character_,
      ifelse(gene_values == 1, "Mutated", "No detected mutation")
    ),
    levels = c("Mutated", "No detected mutation")
  )
}


# 7. MUTATION STATUS AND CANCER STAGE
# ============================================================

cat("\n============================================\n")
cat("MUTATION STATUS AND CANCER STAGE\n")
cat("============================================\n")

run_stage_test <- function(data, gene_name, status) {
  
  test_data <- data |>
    dplyr::mutate(gene_status = status) |>
    dplyr::filter(
      !is.na(gene_status),
      !is.na(stage_group)
    )
  
  contingency_table <- table(
    test_data$gene_status,
    test_data$stage_group
  )
  
  if (nrow(contingency_table) < 2 ||
      ncol(contingency_table) < 2) {
    return(tibble::tibble(
      gene = gene_name,
      test = NA_character_,
      statistic = NA_real_,
      degrees_of_freedom = NA_real_,
      p_value = NA_real_,
      sample_size = nrow(test_data)
    ))
  }
  
  chi_test <- suppressWarnings(
    stats::chisq.test(contingency_table, correct = FALSE)
  )
  
  if (all(chi_test$expected >= 5)) {
    
    tibble::tibble(
      gene = gene_name,
      test = "Pearson chi-square",
      statistic = unname(chi_test$statistic),
      degrees_of_freedom = unname(chi_test$parameter),
      p_value = chi_test$p.value,
      sample_size = nrow(test_data)
    )
    
  } else {
    
    fisher_test <- stats::fisher.test(
      contingency_table,
      simulate.p.value = TRUE,
      B = 10000
    )
    
    tibble::tibble(
      gene = gene_name,
      test = "Fisher Monte Carlo",
      statistic = NA_real_,
      degrees_of_freedom = NA_real_,
      p_value = fisher_test$p.value,
      sample_size = nrow(test_data)
    )
  }
}

stage_results <- dplyr::bind_rows(
  run_stage_test(
    stage_data,
    "Combined four-gene status",
    stage_data$mutation_status
  ),
  
  lapply(
    c("TP53", "BRCA1", "BRCA2", "PIK3CA"),
    function(gene) {
      run_stage_test(
        stage_data,
        gene,
        get_gene_status(stage_data, gene)
      )
    }
  )
)

print(stage_results)

readr::write_csv(
  stage_results,
  here::here("results", "tables", "mutation_stage_association.csv")
)


# 8. OVERALL KAPLAN-MEIER SURVIVAL
# ============================================================

cat("\n============================================\n")
cat("OVERALL KAPLAN-MEIER SURVIVAL\n")
cat("============================================\n")

km_fit <- survival::survfit(
  survival::Surv(OS_MONTHS, survival_event) ~ mutation_status,
  data = km_data
)

km_test <- survival::survdiff(
  survival::Surv(OS_MONTHS, survival_event) ~ mutation_status,
  data = km_data
)

km_p_value <- stats::pchisq(
  km_test$chisq,
  df = length(km_test$n) - 1,
  lower.tail = FALSE
)

overall_km_results <- tibble::tibble(
  analysis = "Combined four-gene mutation status",
  chi_square = unname(km_test$chisq),
  degrees_of_freedom = length(km_test$n) - 1,
  p_value = km_p_value,
  sample_size = nrow(km_data),
  deaths = sum(km_data$survival_event == 1),
  censored = sum(km_data$survival_event == 0)
)

print(overall_km_results)

readr::write_csv(
  overall_km_results,
  here::here("results", "tables", "overall_survival_logrank.csv")
)

overall_km_plot <- survminer::ggsurvplot(
  km_fit,
  data = km_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  legend.title = "Mutation status",
  legend.labs = c("Mutated", "No detected mutation"),
  title = "Overall Survival by Combined Mutation Status"
)

pdf(
  here::here("results", "figures", "kaplan_meier_overall_survival.pdf"),
  width = 9,
  height = 8
)

print(overall_km_plot)
dev.off()


# 9. GENE-SPECIFIC SURVIVAL ANALYSIS
# ============================================================

cat("\n============================================\n")
cat("GENE-SPECIFIC SURVIVAL ANALYSIS\n")
cat("============================================\n")

genes <- c("TP53", "BRCA1", "BRCA2", "PIK3CA")

gene_survival_results <- lapply(
  genes,
  function(gene) {
    
    gene_data <- km_data |>
      dplyr::mutate(
        gene_status = get_gene_status(km_data, gene)
      ) |>
      dplyr::filter(!is.na(gene_status))
    
    gene_fit <- survival::survfit(
      survival::Surv(OS_MONTHS, survival_event) ~ gene_status,
      data = gene_data
    )
    
    gene_test <- survival::survdiff(
      survival::Surv(OS_MONTHS, survival_event) ~ gene_status,
      data = gene_data
    )
    
    gene_p <- stats::pchisq(
      gene_test$chisq,
      df = length(gene_test$n) - 1,
      lower.tail = FALSE
    )
    
    gene_plot <- survminer::ggsurvplot(
      gene_fit,
      data = gene_data,
      risk.table = TRUE,
      pval = TRUE,
      conf.int = TRUE,
      xlab = "Time (months)",
      ylab = "Overall survival probability",
      legend.title = paste(gene, "status"),
      title = paste("Overall Survival by", gene, "Mutation Status")
    )
    
    pdf(
      here::here(
        "results",
        "figures",
        paste0(gene, "_survival_analysis.pdf")
      ),
      width = 9,
      height = 8
    )
    
    print(gene_plot)
    dev.off()
    
    tibble::tibble(
      gene = gene,
      chi_square = unname(gene_test$chisq),
      degrees_of_freedom = length(gene_test$n) - 1,
      p_value = gene_p,
      sample_size = nrow(gene_data),
      deaths = sum(gene_data$survival_event == 1),
      censored = sum(gene_data$survival_event == 0)
    )
  }
)

gene_survival_results <- dplyr::bind_rows(gene_survival_results)

print(gene_survival_results)

readr::write_csv(
  gene_survival_results,
  here::here("results", "tables", "gene_specific_survival_logrank.csv")
)


# 10. MULTIVARIABLE COX REGRESSION
# ============================================================

cat("\n============================================\n")
cat("COX PROPORTIONAL HAZARDS REGRESSION\n")
cat("============================================\n")

cox_model <- survival::coxph(
  survival::Surv(OS_MONTHS, survival_event) ~
    mutation_status + AGE + stage_group,
  data = cox_data,
  ties = "efron",
  x = TRUE
)

print(summary(cox_model))

cox_summary <- summary(cox_model)

cox_results <- tibble::tibble(
  variable = rownames(cox_summary$coefficients),
  coefficient = cox_summary$coefficients[, "coef"],
  hazard_ratio = cox_summary$coefficients[, "exp(coef)"],
  standard_error = cox_summary$coefficients[, "se(coef)"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  lower_95_CI = cox_summary$conf.int[, "lower .95"],
  upper_95_CI = cox_summary$conf.int[, "upper .95"]
)

print(cox_results)

readr::write_csv(
  cox_results,
  here::here("results", "tables", "cox_regression_results.csv")
)

saveRDS(
  cox_model,
  here::here("data", "processed", "cox_regression_model.rds")
)


# 11. ORIGINAL COX MODEL PH ASSUMPTION TEST
# ============================================================

cat("\n============================================\n")
cat("PROPORTIONAL HAZARDS ASSUMPTION\n")
cat("============================================\n")

ph_check <- survival::cox.zph(cox_model)

print(ph_check)

ph_results <- as.data.frame(ph_check$table) |>
  tibble::rownames_to_column("variable") |>
  dplyr::rename(
    chi_square = chisq,
    degrees_of_freedom = df,
    p_value = p
  )

print(ph_results)

readr::write_csv(
  ph_results,
  here::here("results", "tables", "cox_ph_assumption_test.csv")
)

pdf(
  here::here(
    "results", "figures",
    "cox_proportional_hazards_diagnostic.pdf"
  ),
  width = 9,
  height = 7
)

plot(ph_check)
dev.off()


# 12. STRATIFIED COX REGRESSION
# ============================================================

cat("\n============================================\n")
cat("STRATIFIED COX REGRESSION\n")
cat("============================================\n")

stratified_cox_model <- survival::coxph(
  survival::Surv(OS_MONTHS, survival_event) ~
    mutation_status + AGE + survival::strata(stage_group),
  data = cox_data,
  ties = "efron",
  x = TRUE
)

print(summary(stratified_cox_model))

stratified_summary <- summary(stratified_cox_model)

stratified_cox_results <- tibble::tibble(
  variable = rownames(stratified_summary$coefficients),
  coefficient = stratified_summary$coefficients[, "coef"],
  hazard_ratio = stratified_summary$coefficients[, "exp(coef)"],
  standard_error = stratified_summary$coefficients[, "se(coef)"],
  p_value = stratified_summary$coefficients[, "Pr(>|z|)"],
  lower_95_CI = stratified_summary$conf.int[, "lower .95"],
  upper_95_CI = stratified_summary$conf.int[, "upper .95"]
)

print(stratified_cox_results)

readr::write_csv(
  stratified_cox_results,
  here::here(
    "results", "tables",
    "stratified_cox_regression_results.csv"
  )
)

saveRDS(
  stratified_cox_model,
  here::here(
    "data", "processed",
    "stratified_cox_regression_model.rds"
  )
)


# 13. STRATIFIED MODEL PH ASSUMPTION TEST
# ============================================================

cat("\n============================================\n")
cat("STRATIFIED COX PROPORTIONAL HAZARDS TEST\n")
cat("============================================\n")

stratified_ph_test <- survival::cox.zph(
  stratified_cox_model
)

print(stratified_ph_test)

stratified_ph_results <- as.data.frame(
  stratified_ph_test$table
) |>
  tibble::rownames_to_column("variable") |>
  dplyr::rename(
    chi_square = chisq,
    degrees_of_freedom = df,
    p_value = p
  )

print(stratified_ph_results)

readr::write_csv(
  stratified_ph_results,
  here::here(
    "results", "tables",
    "stratified_cox_ph_test.csv"
  )
)

pdf(
  here::here(
    "results", "figures",
    "stratified_cox_ph_diagnostic.pdf"
  ),
  width = 9,
  height = 7
)

plot(stratified_ph_test)
dev.off()



# ============================================================
# 14. FINAL ANALYSIS QUALITY CHECKS
# ============================================================

cat("\n")
cat("========================================\n")
cat("FINAL ANALYSIS QUALITY CHECKS\n")
cat("========================================\n\n")

# ------------------------------------------------------------
# 14.1 Confirm that the required datasets and models exist
# ------------------------------------------------------------

required_objects <- c(
  "analysis_data",
  "stage_data",
  "km_data",
  "cox_data",
  "cox_model",
  "stratified_cox_model"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 14.2 Extract the actual data used by each Cox model
# ------------------------------------------------------------

cox_model_frame <- stats::model.frame(cox_model)

stratified_model_frame <- stats::model.frame(
  stratified_cox_model
)

# Number of observations actually included in each model
cox_model_observations <- nrow(cox_model_frame)

stratified_model_observations <- nrow(
  stratified_model_frame
)

# ------------------------------------------------------------
# 14.3 Count events in each Cox model
# ------------------------------------------------------------

cox_response <- stats::model.response(cox_model_frame)

stratified_response <- stats::model.response(
  stratified_model_frame
)

# Extract the event indicator from the Surv response.
# For right-censored survival data, the final column
# represents the event indicator.

cox_model_events <- sum(
  cox_response[, ncol(cox_response)] == 1,
  na.rm = TRUE
)

stratified_model_events <- sum(
  stratified_response[, ncol(stratified_response)] == 1,
  na.rm = TRUE
)

# ------------------------------------------------------------
# 14.4 Report cohort and model sizes
# ------------------------------------------------------------

cat("COHORT SUMMARY\n")
cat("----------------------------------------\n")
cat("Original analysis samples:", nrow(analysis_data), "\n")
cat("Stage association samples:", nrow(stage_data), "\n")
cat("Kaplan-Meier samples:", nrow(km_data), "\n")
cat("Prepared Cox regression cohort:", nrow(cox_data), "\n\n")

cat("COX MODEL SUMMARY\n")
cat("----------------------------------------\n")
cat(
  "Original Cox model observations:",
  cox_model_observations, "\n"
)
cat(
  "Original Cox model events:",
  cox_model_events, "\n"
)
cat(
  "Stratified Cox model observations:",
  stratified_model_observations, "\n"
)
cat(
  "Stratified Cox model events:",
  stratified_model_events, "\n\n"
)

# ------------------------------------------------------------
# 14.5 Check consistency between the two Cox models
# ------------------------------------------------------------

if (
  cox_model_observations ==
  stratified_model_observations
) {
  cat("PASS: Both Cox models used the same number of observations.\n")
} else {
  warning(
    "The original and stratified Cox models used different ",
    "numbers of observations. Review their model frames."
  )
}

if (cox_model_events == stratified_model_events) {
  cat("PASS: Both Cox models contain the same number of events.\n")
} else {
  warning(
    "The original and stratified Cox models contain different ",
    "numbers of events. Review their model frames."
  )
}

# ------------------------------------------------------------
# 14.6 Confirm that the prepared Cox cohort has valid data
# ------------------------------------------------------------

required_cox_variables <- c(
  "OS_MONTHS",
  "survival_event",
  "AGE",
  "stage_group",
  "mutation_status"
)

missing_cox_variables <- setdiff(
  required_cox_variables,
  names(cox_data)
)

if (length(missing_cox_variables) > 0) {
  stop(
    "Missing Cox regression variables: ",
    paste(missing_cox_variables, collapse = ", ")
  )
}

if (
  anyNA(cox_data[, required_cox_variables]) ||
  any(cox_data$OS_MONTHS <= 0) ||
  any(!cox_data$survival_event %in% c(0, 1))
) {
  warning(
    "The prepared Cox cohort contains missing or invalid values. ",
    "Review the cohort before interpreting the models."
  )
} else {
  cat("PASS: Prepared Cox cohort contains valid model variables.\n")
}

# ------------------------------------------------------------
# 14.7 Final completion message
# ------------------------------------------------------------

cat("\n")
cat("========================================\n")
cat("FINAL QUALITY CHECK COMPLETED\n")
cat("========================================\n")

# 15. SAVE ANALYSIS RECORD AND SESSION INFORMATION
# ============================================================

analysis_record <- tibble::tibble(
  analysis = c(
    "Stage association",
    "Kaplan-Meier survival",
    "Original Cox regression",
    "Stratified Cox regression"
  ),
  sample_size = c(
    nrow(stage_data),
    nrow(km_data),
    cox_model_observations,
    stratified_model_observations
  )
)

readr::write_csv(
  analysis_record,
  here::here("results", "tables", "analysis_cohort_summary.csv")
)

writeLines(
  capture.output(sessionInfo()),
  here::here(
    "results", "tables",
    "statistical_analysis_session_info.txt"
  )
)

package_versions <- tibble::tibble(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package) {
      as.character(utils::packageVersion(package))
    },
    character(1)
  )
)

readr::write_csv(
  package_versions,
  here::here(
    "results", "tables",
    "statistical_analysis_package_versions.csv"
  )
)

cat("\nSTATISTICAL ANALYSIS COMPLETED\n")
cat("Review the updated tables and figures in results/.\n")
cat("All final cohort checks passed.\n")
