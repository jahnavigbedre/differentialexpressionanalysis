# ============================================================
# Script: edgeR_multiMiR_pipeline.R
# Purpose: Differential expression analysis with edgeR 
#          followed by validated target retrieval using multiMiR.
# ============================================================

# ---- Load Required Packages ----
packages <- c("edgeR", "multiMiR")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ---- Input Data ----
# Example files provided within multiMiR package
counts_file  <- system.file("extdata", "counts_table.tsv", package = "multiMiR")
strains_file <- system.file("extdata", "strains_factor.tsv", package = "multiMiR")

# Read input files
counts_table   <- readRDS(counts_file)
strains_factor <- readRDS(strains_file)

# ---- Experimental Design ----
design <- model.matrix(~ strains_factor)

# ---- edgeR Differential Expression Analysis ----
dge <- DGEList(counts = counts_table)
dge <- calcNormFactors(dge)
dge$samples$strains <- strains_factor

# Estimate dispersions
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# Fit GLM model
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Extract DE results
p_val_DE_edgeR <- topTags(lrt, adjust.method = "BH", n = Inf)

# ---- Top Differentially Expressed miRNAs ----
top_miRNAs <- rownames(p_val_DE_edgeR$table)[1:10]
cat("✅ Top 10 differentially expressed miRNAs:\n")
print(top_miRNAs)

# ---- Query multiMiR for Validated Targets ----
multimir_results <- get_multimir(
  org     = "mmu",        # mouse example
  mirna   = top_miRNAs,
  table   = "validated",
  summary = TRUE
)

# Display results
cat("✅ multiMiR results retrieved:\n")
print(head(multimir_results@data))
f
