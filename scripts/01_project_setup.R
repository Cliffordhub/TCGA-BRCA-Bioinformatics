
# ============================================================
# PROJECT: TCGA-BRCA Bioinformatics
# SCRIPT 01: Project Setup and Environment Preparation
# ============================================================

# 1. Confirm the project working directory -------------------

project_root <- getwd()

cat("\n============================================\n")
cat("TCGA-BRCA PROJECT SETUP\n")
cat("============================================\n")

cat("Project directory:\n", project_root, "\n")
cat("R version:", R.version.string, "\n\n")

# Check that this is the correct project directory
if (!dir.exists(file.path(project_root, "scripts"))) {
  stop(
    "The scripts folder was not found. ",
    "Open the TCGA-BRCA-Bioinformatics RStudio project ",
    "before running this script."
  )
}

# 2. Create the project folder structure ---------------------

project_folders <- c(
  "data/raw",
  "data/processed",
  "scripts",
  "results/figures",
  "results/tables",
  "report"
)

for (folder in project_folders) {
  dir.create(
    file.path(project_root, folder),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("Project folders checked/created.\n")

# 3. Check and install required packages ---------------------

cran_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "survival",
  "survminer",
  "here",
  "janitor",
  "tibble"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace,
          logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  cat("\nInstalling missing CRAN packages:\n")
  print(missing_cran)
  
  install.packages(missing_cran)
} else {
  cat("All required CRAN packages are installed.\n")
}

# cBioPortalData is distributed through Bioconductor
if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
  
  cat("\nChecking BiocManager installation...\n")
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  cat("Installing cBioPortalData...\n")
  
  BiocManager::install(
    "cBioPortalData",
    ask = FALSE,
    update = FALSE
  )
  
} else {
  cat("cBioPortalData is already installed.\n")
}

# 4. Load the packages used in the project -------------------

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "survival",
  "survminer",
  "here",
  "janitor",
  "tibble",
  "cBioPortalData"
)

for (package in required_packages) {
  suppressPackageStartupMessages(
    library(package, character.only = TRUE)
  )
}

cat("\nRequired packages loaded successfully.\n")

# 5. Record package versions ---------------------------------

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

cat("\nPACKAGE VERSION CHECK\n")
print(package_versions, row.names = FALSE)

# Save package versions for project documentation
write.csv(
  package_versions,
  file = file.path(
    project_root,
    "results/tables/package_versions.csv"
  ),
  row.names = FALSE
)

# 6. Record the R session information ------------------------

writeLines(
  capture.output(sessionInfo()),
  con = file.path(
    project_root,
    "results/tables/session_info_setup.txt"
  )
)

# 7. Final setup quality checks ------------------------------

stopifnot(
  all(dir.exists(file.path(project_root, project_folders))),
  all(vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  ))
)

cat("\n============================================\n")
cat("SCRIPT 01 COMPLETED SUCCESSFULLY\n")
cat("Project environment is ready.\n")
cat("============================================\n")
