# ============================================================
# Script: repeated_genes_finder.R
# Purpose: Identify unique string values that are repeated 
#          across all (or a threshold of) columns in an Excel file.
# ============================================================

# ---- Install & Load Packages ----
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

# ---- Parameters (User-Defined) ----
input_file <- "Copy of CDX_only.xlsx"   # Path to Excel file
sheet_name <- 1                         # Sheet index or name
threshold <- 9                          # Minimum number of occurrences across dataset
output_file <- "Repeated_Genes.csv"     # Output file (CSV)

# ---- Read Data ----
df <- read_excel(input_file, sheet = sheet_name)

# ---- Function: Check if value is repeated ----
is_repeated_in_all_columns <- function(value, threshold) {
  all_values <- unlist(df)
  num_occurrences <- sum(all_values == value, na.rm = TRUE)
  num_occurrences >= threshold   # Condition: Appears at least 'threshold' times
}

# ---- Find Unique Values ----
unique_values <- unique(unlist(df))
repeated_strings <- unique_values[sapply(unique_values, is_repeated_in_all_columns, threshold = threshold)]

# ---- Sort Results ----
sorted_repeated_strings <- sort(repeated_strings)

# ---- Convert to DataFrame ----
result_df <- data.frame(Sorted_Repeated_Genes = sorted_repeated_strings)

# ---- Save Results to CSV ----
write.csv(result_df, output_file, row.names = FALSE)

# ---- Print Output ----
cat("✅ Repeated values saved to:", output_file, "\n")
print(result_df)
fyour
