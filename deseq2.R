# ============================================================
# Script: deseq2_analysis_pipeline.R
# Purpose: Differential expression analysis using DESeq2
# ============================================================

# ---- Load Required Packages ----
packages <- c("DESeq2", "RColorBrewer", "gplots", "pheatmap", "genefilter", "calibrate")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("DESeq2")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ---- Input Data ----
countdata <- read.table("counts.txt", header = TRUE, row.names = 1)

# Keep relevant columns (remove first 5 if they are annotation)
countdata <- countdata[, 6:ncol(countdata)]

# Clean column names
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix
countdata <- as.matrix(countdata)
head(countdata)

# ---- Experimental Conditions ----
condition <- factor(c(rep("ctl", 74), rep("dis", 74)))
coldata <- data.frame(row.names = colnames(countdata), condition)

# ---- DESeq2 Analysis ----
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~condition)
dds <- DESeq(dds)

# ---- QC Plots ----
png("qc-dispersions.png", 1000, 1000, pointsize = 20)
plotDispEsts(dds, main = "Dispersion plot")
dev.off()

# Regularized log transformation
rld <- rlogTransformation(dds)

# Heatmap of sample distances
mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
sampleDists <- as.matrix(dist(t(assay(rld))))

png("qc-heatmap-samples.png", w = 1000, h = 1000, pointsize = 20)
heatmap.2(as.matrix(sampleDists), key = FALSE, trace = "none",
          col = colorpanel(100, "black", "white"),
          ColSideColors = mycols[condition], RowSideColors = mycols[condition],
          margin = c(10, 10), main = "Sample Distance Matrix")
dev.off()

# Top 100 expressed genes heatmap
select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:100]
df <- as.data.frame(colData(dds)[, c("condition", "sizeFactor")])
ntd <- normTransform(dds)

pheatmap(assay(ntd)[select, ],
         cluster_rows = TRUE, show_rownames = TRUE,
         cluster_cols = FALSE, annotation_col = df,
         cutree_rows = 1,
         color = colorRampPalette(c("blue", "white", "red"))(50))

# ---- PCA Plot ----
rld_pca <- function(rld, intgroup = "condition", ntop = 500,
                    colors = NULL, legendpos = "bottomleft",
                    main = "PCA Biplot", textcx = 1, ...) {
  rv <- rowVars(assay(rld))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(rld)[select, ]))
  fac <- factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))

  if (is.null(colors)) {
    colors <- if (nlevels(fac) >= 3) brewer.pal(nlevels(fac), "Paired") else c("black", "red")
  }

  pc1var <- round(summary(pca)$importance[2, 1] * 100, digits = 1)
  pc2var <- round(summary(pca)$importance[2, 2] * 100, digits = 1)
  pc1lab <- paste0("PC1 (", pc1var, "%)")
  pc2lab <- paste0("PC2 (", pc2var, "%)")

  plot(PC2 ~ PC1, data = as.data.frame(pca$x), bg = colors[fac],
       pch = 21, xlab = pc1lab, ylab = pc2lab, main = main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs = rownames(as.data.frame(pca$x)), cex = textcx))
  legend(legendpos, legend = levels(fac), col = colors, pch = 20)
}

png("qc-pca.png", 1000, 1000, pointsize = 20)
rld_pca(rld, colors = mycols, intgroup = "condition", xlim = c(-75, 35))
dev.off()

# ---- Differential Expression Results ----
res <- results(dds)
table(res$padj < 0.05)

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Merge with normalized counts
resdata <- merge(as.data.frame(res),
                 as.data.frame(counts(dds, normalized = TRUE)),
                 by = "row.names", sort = FALSE)
names(resdata)[1] <- "Gene"

# Save results
write.csv(resdata, file = "diffexpr-results.csv")

# ---- Diagnostic Plots ----
hist(res$pvalue, breaks = 50, col = "grey", main = "P-value distribution")

attr(res, "filterThreshold")
plot(attr(res, "filterNumRej"), type = "b",
     xlab = "Quantiles of baseMean", ylab = "Number of rejections")

# ---- MA Plot ----
maplot <- function(res, thresh = 0.05, labelsig = FALSE, textcx = 1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch = 20, cex = .5, log = "x", ...))
  with(subset(res, padj < thresh), points(baseMean, log2FoldChange, col = "red", pch = 20, cex = 1.5))
  if (labelsig) {
    with(subset(res, padj < thresh), textxy(baseMean, log2FoldChange, labs = Gene, cex = textcx, col = 2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize = 20)
maplot(resdata, main = "MA Plot")
dev.off()

# ---- Volcano Plot ----
volcanoplot <- function(res, lfcthresh = 2, sigthresh = 0.05,
                        main = "Volcano Plot", legendpos = "bottomright",
                        labelsig = TRUE, textcx = 1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch = 20, main = main, ...))
  with(subset(res, padj < sigthresh), points(log2FoldChange, -log10(pvalue), col = "red", pch = 20))
  with(subset(res, abs(log2FoldChange) > lfcthresh), points(log2FoldChange, -log10(pvalue), col = "orange", pch = 20))
  with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh),
       points(log2FoldChange, -log10(pvalue), col = "green", pch = 20))
  if (labelsig) {
    with(subset(res, padj < sigthresh & abs(log2FoldChange) > lfcthresh),
         textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = textcx))
  }
  legend(legendpos, legend = c(paste("FDR <", sigthresh), paste("|LogFC| >", lfcthresh), "Both"),
         pch = 20, col = c("red", "orange", "green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize = 20)
volcanoplot(resdata, lfcthresh = 1, sigthresh = 0.05, textcx = .8, xlim = c(-2.3, 2))
dev.off()

cat("✅ DESeq2 pipeline completed. Results saved to diffexpr-results.csv\n")
