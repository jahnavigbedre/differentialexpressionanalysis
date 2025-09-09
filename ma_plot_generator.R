# ============================================================
# Script: volcano_plot_ncRNA.R
# Purpose: Generate volcano plot from Excel-based differential 
#          expression data using ggmaplot.
# ============================================================

# ---- Install & Load Required Packages ----
packages <- c("readxl", "ggplot2", "reshape2", "ggpubr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ---- Parameters ----
input_file <- "E:/Desktop 23/Articles in Processing/ncRNA/REVISED/Book1.xlsx"
sheet_name <- "Sheet2"
gene_column <- "name"   # Column with gene names
fdr_threshold <- 0.05
fc_threshold <- 2

# Custom colors for up, down, and non-significant
my_palette <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "gray70")

# ---- Read Data ----
data <- read_excel(input_file, sheet = sheet_name)

# ---- Check data ----
head(data)

# ---- Volcano Plot ----
ggmaplot(
  data,
  main = expression("Group 1" %->% "Group 2"),
  fdr = fdr_threshold,
  fc = fc_threshold,
  size = 0.6,
  palette = my_palette,
  genenames = as.vector(data[[gene_column]]),
  legend = "top",
  top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = ggplot2::theme_minimal()
)
