# 🧬 miRNA-DE-multiMiR

<p align="center">
  <img src="https://img.shields.io/badge/R-4.3%2B-blue?style=flat-square&logo=r" />
  <img src="https://img.shields.io/badge/edgeR-DE%20Analysis-green?style=flat-square" />
  <img src="https://img.shields.io/badge/multiMiR-Target%20Validation-purple?style=flat-square" />
  <img src="https://img.shields.io/badge/License-MIT-yellow?style=flat-square" />
</p>

---
# edgeR + multiMiR Pipeline

This repository contains an R workflow for performing **differential expression analysis** of small RNAs (e.g., miRNAs) using **edgeR**, followed by retrieval of **validated miRNA targets** with the [multiMiR](https://cran.r-project.org/web/packages/multiMiR/) package.

---

## 📌 Features
- Reads count data and experimental conditions.  
- Performs **differential expression analysis** using `edgeR`.  
- Extracts top differentially expressed miRNAs.  
- Retrieves **validated miRNA targets** from the multiMiR database.  
- Outputs results for downstream biological interpretation.  

---

## 📦 Requirements

Make sure you have the following R packages installed:

- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)  
- [multiMiR](https://cran.r-project.org/web/packages/multiMiR/)  

To install them:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor package
BiocManager::install("edgeR")

# CRAN package
install.packages("multiMiR")
````

---

## 📂 Input Files

The pipeline requires:

1. **Counts table** (`counts_table.tsv`) – matrix of raw read counts (genes/miRNAs as rows, samples as columns).
2. **Strains/condition factor** (`strains_factor.tsv`) – metadata describing sample groups.

👉 Example datasets are included with the `multiMiR` package (`extdata/`).

---

## 🚀 Usage

Run the script:

```r
source("edgeR_multiMiR_pipeline.R")
```

This will:

1. Load the count data and experimental design.
2. Perform normalization, dispersion estimation, and differential expression analysis with edgeR.
3. Extract the **top 10 differentially expressed miRNAs**.
4. Query multiMiR for validated targets in mouse (`org = "mmu"`).

---

## 📊 Output

* **Top 10 DE miRNAs** are printed in the console.
* **multiMiR target results** are displayed (and can be saved).

Example console output:

```
✅ Top 10 differentially expressed miRNAs:
[1] "mmu-miR-21a-5p" "mmu-miR-155-5p" "mmu-miR-34a-5p" ...

✅ multiMiR results retrieved:
       mature_mirna_id target_symbol database
1       mmu-miR-21a-5p          PTEN    miRTarBase
2       mmu-miR-155-5p           SOCS1  TarBase
...
```

---

## 🧬 Customization

* Change `org = "mmu"` to another organism (e.g., `"hsa"` for human).
* Adjust the number of top DE miRNAs (`[1:10]`) as needed.
* Save results with:

```r
write.csv(p_val_DE_edgeR$table, "DE_miRNAs_edgeR.csv")
write.csv(multimir_results@data, "multiMiR_targets.csv")
```

---

## 📖 References

* **edgeR**: Robinson MD, McCarthy DJ, Smyth GK (2010). *edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.* Bioinformatics, 26(1), 139–140.
* **multiMiR**: Ru Y, Kechris KJ, Tabakoff B, et al. (2014). *The multiMiR R package and database: integration of microRNA–target interactions along with their disease and drug associations.* Nucleic Acids Research, 42(17).

