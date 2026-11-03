
# ==========================================
# 01. INSTALLATION OF ANALYSIS PACKAGES
# ==========================================

# Packages needed for our research
packages <- c(
  "tidyverse",
  "survminer",
  "janitor",
  "here",
  "rmarkdown"
)

# Check which packages are missing
new_packages <- packages[
  !packages %in% installed.packages()[, "Package"]
]

# Install only missing packages
if (length(new_packages) > 0) {
  install.packages(new_packages)
}

# Load the packages
library(tidyverse)
library(survminer)
library(janitor)
library(here)
library(rmarkdown)

# Check installation
print("Analysis packages loaded successfully!")
sessionInfo()
R.version.string
