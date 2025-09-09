# ============================================================
# Script: install_required_packages.R
# Purpose: Install all required CRAN and Bioconductor packages 
#          for differential expression and visualization in R.
# ============================================================

# ---- Install BiocManager (if not already installed) ----
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ---- CRAN Packages ----
cran_packages <- c("ggplot2", "RColorBrewer", "gplots", "genefilter")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# ---- Bioconductor Packages ----
bioc_packages <- c("DESeq2")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# ---- Load All Packages ----
lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)

cat("✅ All required packages are installed and loaded!\n")
