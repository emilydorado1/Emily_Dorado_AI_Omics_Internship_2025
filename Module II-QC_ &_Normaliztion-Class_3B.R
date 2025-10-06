#-------------
# Assignment- Module II-QC_&_Normalization_Class_3B
# Emily Dorado
#-------------

# Work through the full preprocessing workflow with your own dataset
# Dataset I'm working with: GSE108134
# Focuses 

# First install all nessecary/required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("Biobase","limma","arrayQualityMetrics","matrixStats")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p, ask=FALSE)

cran_pkgs <- c("dplyr","readr","stringr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

install.packages(c("data.table", "stringr", "matrixStats", "readr", "dplyr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("arrayQualityMetrics", ask = FALSE, update = TRUE) 

library(Biobase)
library(limma)
library(matrixStats)
library(dplyr)
library(readr)
library(stringr)


set.seed(108134)

# Load in my database
# File: GSE108134_307SAE-S_PM2.5_resid.csv  (tab-delimited)
# install.packages(c("readr","data.table"), repos="https://cloud.r-project.org")
library(readr); library(data.table)

path <- file.choose()  # pick your file in the dialog
# Try TSV first (your file is TSV). If that fails, try CSV.
expr_df <- tryCatch(
  read_tsv(path, show_col_types = FALSE),
  error = function(e) tryCatch(read_csv(path, show_col_types = FALSE), error = function(e2) fread(path) |> as.data.frame())
)
colnames(expr_df)[1] <- "ID"  # ensure first column is named ID

# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization

# Safety checks
stopifnot(colnames(expr_df)[1] == "ID")
expr_df <- expr_df[!duplicated(expr_df$ID), ]     # ensure unique probe IDs

# ---- Build ExpressionSet ------------------------------------------------------
library(Biobase); library(stringr); library(matrixStats)

probe_ids <- expr_df$ID
expr_mat  <- as.matrix(expr_df[,-1])
mode(expr_mat) <- "numeric"
samples <- colnames(expr_mat)

# Working on Pheno Type
# Expected sample name pattern: DGM-12345_M0_sm  (subject, timepoint, group)
subject   <- stringr::str_extract(samples, "^DGM-[0-9]+")
timepoint <- stringr::str_extract(samples, "_M[0-9]+_") |> stringr::str_replace_all("[_M]", "")
group     <- stringr::str_extract(samples, "_[A-Za-z]+$") |> stringr::str_replace_all("_", "")

pheno <- data.frame(
  Sample = samples,
  Subject = subject,
  Timepoint = factor(timepoint, levels = c("0","3","6","12"), labels = c("M0","M3","M6","M12")),
  Group = factor(group),  # in this file: likely all "sm"
  BaselineVsFollowup = factor(ifelse(timepoint == "0", "baseline", "followup"),
                              levels = c("baseline","followup")),
  row.names = samples,
  check.names = FALSE
)

fData_df <- data.frame(ID = probe_ids, row.names = probe_ids)

eset_raw <- ExpressionSet(
  assayData  = expr_mat,
  phenoData  = AnnotatedDataFrame(pheno),
  featureData= AnnotatedDataFrame(fData_df)
)

cat("Samples:", ncol(exprs(eset_raw)), " | Probes:", nrow(exprs(eset_raw)), "\n")
cat("Timepoints:", paste(levels(pheno$Timepoint), table(pheno$Timepoint), sep=":", collapse=" | "), "\n")
cat("Group counts:", paste(levels(pheno$Group), table(pheno$Group), sep=":", collapse=" | "), "\n")

#----------------------------------------------------------------------
# 1) QC BEFORE normalization
# ---------------------------------------------------------------------
# Create directory to store plots
dir.create("Results/QC_Pre", recursive = TRUE, showWarnings = FALSE)

# Extract expression matrix
X <- exprs(eset_raw)

# ----------------------------------------------------------
#### Boxplot: Raw intensity distribution across samples ####
# ----------------------------------------------------------
pdf("Results/QC_Pre/Boxplot_Raw_Data.pdf", width = 12, height = 6)
par(mar = c(7,4,2,1))
boxplot(as.data.frame(X),
        las = 2, outline = FALSE,
        main = "Pre-Normalization: Boxplot per Array",
        ylab = "Expression Intensity")
abline(h = median(X, na.rm = TRUE), lty = 2, col = "gray40")
dev.off()

# ----------------------------------------------------------
#### Density Plot: Visualize overall intensity distribution ####
# ----------------------------------------------------------
pdf("Results/QC_Pre/Density_Raw_Data.pdf", width = 8, height = 6)
limma::plotDensities(X,
                     legend = FALSE,
                     main = "Pre-Normalization: Density Plot")

if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Results/QC_Post")) dir.create("Results/QC_Post")

# Requires: limma
if (!requireNamespace("limma", quietly = TRUE)) install.packages("limma")
library(limma)

# eset_raw must already exist (built from your expr_df)
stopifnot(exists("eset_raw"))

# Quantile normalize the expression matrix from eset_raw
X_raw <- exprs(eset_raw)
X_norm <- normalizeBetweenArrays(X_raw, method = "quantile")

# Create normalized ExpressionSet
eset_norm <- eset_raw
exprs(eset_norm) <- X_norm


pdf("Results/QC_Post/Density_Normalized_Data.pdf", width = 8, height = 6)
limma::plotDensities(X_norm,
                     legend = FALSE,
                     main = "Post-Normalization: Density Plot")

# --- PDF boxplot ---
pdf("Results/QC_Post/Boxplot_Normalized_Data.pdf", width = 12, height = 6)
par(mar = c(7,4,2,1))
boxplot(as.data.frame(X_norm),
        las = 2, outline = FALSE,
        main = "Post-Normalization: Boxplot per Array",
        ylab = "Normalized Expression")
abline(h = median(X_norm, na.rm = TRUE), lty = 2, col = "gray40")
dev.off()

# --- PNG boxplot (easy to upload to the form) ---
png("Results/QC_Post/Boxplot_Normalized_Data.png", width = 1600, height = 900, res = 150)
par(mar = c(7,4,2,1))
boxplot(as.data.frame(X_norm),
        las = 2, outline = FALSE,
        main = "Post-Normalization: Boxplot per Array",
        ylab = "Normalized Expression")
abline(h = median(X_norm, na.rm = TRUE), lty = 2, col = "gray40")
dev.off()

#-----------------------------------------------------------
# Correlation Heatmap: Check sample similarity 
# ----------------------------------------------------------
cor_mat <- cor(X, use = "pairwise.complete.obs")
pdf("Results/QC_Pre/Correlation_Heatmap_Raw_Data.pdf", width = 7, height = 7)
heatmap(cor_mat, symm = TRUE,
        main = "Pre-Normalization: Sample–Sample Correlation")
dev.off()

# ----------------------------------------------------------
# Outlier Flagging (simple heuristic)
# ----------------------------------------------------------
# Identify samples with low median correlation to others
med_corr <- apply(cor_mat, 2,
                  function(v) median(v[!is.na(v) & v < 0.9999]))
corr_threshold <- median(med_corr) - 3 * mad(med_corr)
corr_outlier_flag <- med_corr < corr_threshold | med_corr < 0.9

# Store QC results
qc_summary <- data.frame(
  Sample = names(med_corr),
  MedianCorrelation = med_corr,
  Outlier = corr_outlier_flag
)
write.csv(qc_summary,
          "Results/QC_Pre/QC_Raw_Summary.csv",
          row.names = FALSE)
# Count number of flagged arrays
pre_outliers <- sum(qc_summary$Outlier)
cat("Pre-Normalization outliers flagged:", pre_outliers, "\n")

#### Note: I made the boxplot and it look horrible with all the data in it 
### Goal: Make it look much nicer

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
library(ggplot2); library(scales)

outdir <- "Results/QC_Post"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Helper: shorten labels like "DGM-10132_M0_sm" → "10132_M0"
short_labels <- function(x) {
  sub("^DGM-", "", sub("_sm$", "", x))
}

# -------------------------------
# 1) Representative Boxplot (48)
# -------------------------------
set.seed(1)
X <- X_norm
n <- ncol(X)
# choose ~48 evenly spaced indices across all samples
k <- 48
idx <- unique(round(seq(1, n, length.out = k)))
X_rep <- X[, idx, drop = FALSE]

df_rep <- data.frame(
  value = as.vector(X_rep),
  sample = rep(short_labels(colnames(X_rep)), each = nrow(X_rep)),
  stringsAsFactors = FALSE
)

p_box_rep <- ggplot(df_rep, aes(x = sample, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.3) +
  coord_cartesian(ylim = quantile(df_rep$value, c(0.01, 0.99), na.rm = TRUE)) +
  labs(title = "Post-Normalization: Representative Boxplots (48 of 307 arrays)",
       x = "Sample (subset)", y = "Normalized Expression") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(outdir, "Boxplot_Normalized_Representative48.png"),
       p_box_rep, width = 12, height = 5, dpi = 200)
ggsave(file.path(outdir, "Boxplot_Normalized_Representative48.pdf"),
       p_box_rep, width = 12, height = 5)

# -----------------------------------
# 2) Median & IQR per Sample (clean)
# -----------------------------------
med <- apply(X, 2, median, na.rm = TRUE)
q1  <- apply(X, 2, function(v) quantile(v, 0.25, na.rm = TRUE))
q3  <- apply(X, 2, function(v) quantile(v, 0.75, na.rm = TRUE))

df_stats <- data.frame(
  sample = factor(seq_len(n), levels = seq_len(n)),
  label  = short_labels(colnames(X)),
  median = med,
  q1 = q1,
  q3 = q3
)

p_median_iqr <- ggplot(df_stats, aes(x = as.integer(sample))) +
  geom_ribbon(aes(ymin = q1, ymax = q3), alpha = 0.2) +
  geom_line(aes(y = median), linewidth = 0.4) +
  labs(title = "Post-Normalization: Median and IQR per Sample",
       x = "Sample index (1…307)", y = "Normalized Expression") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "Median_IQR_Per_Sample.png"),
       p_median_iqr, width = 12, height = 4, dpi = 200)
ggsave(file.path(outdir, "Median_IQR_Per_Sample.pdf"),
       p_median_iqr, width = 12, height = 4)

# ------------------------------------------------
# 3) Single Summary Boxplot (pooled distribution)
# ------------------------------------------------
df_pooled <- data.frame(value = as.vector(X))
p_box_pooled <- ggplot(df_pooled, aes(x = "All Samples", y = value)) +
  geom_boxplot(width = 0.3, outlier.alpha = 0.2) +
  labs(title = "Post-Normalization: Summary Boxplot (Pooled)",
       x = NULL, y = "Normalized Expression") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 11),
        plot.title = element_text(face = "bold"))

ggsave(file.path(outdir, "Boxplot_Normalized_SummaryPooled.png"),
       p_box_pooled, width = 4, height = 5, dpi = 200)
ggsave(file.path(outdir, "Boxplot_Normalized_SummaryPooled.pdf"),
       p_box_pooled, width = 4, height = 5)

# =====================================================================
# PCA PLOT (After Normalization)
# =====================================================================
# Goal:
#   Visualize clustering or separation among samples after normalization.
# Output:
#   Results/QC_Post/PCA_Normalized_Data.[pdf/png]
# =====================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# PCA on samples (each column = one array)
X <- t(scale(t(X_norm), center = TRUE, scale = FALSE))   # center genes
pca <- prcomp(t(X), scale. = FALSE)
scores <- as.data.frame(pca$x)

# Extract top 2 PCs
var_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
colnames(scores)[1:2] <- c("PC1", "PC2")

# Optional: parse timepoint info from sample names for color labels
samples <- colnames(X_norm)
timepoint <- gsub(".*_M([0-9]+)_.*", "\\1", samples)
timepoint <- factor(timepoint, levels = c("0","3","6","12"),
                    labels = c("M0","M3","M6","M12"))

scores$Timepoint <- timepoint
scores$Sample <- samples

# --- Pretty PCA plot ---
p_pca <- ggplot(scores, aes(x = PC1, y = PC2, color = Timepoint)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("#0072B2","#56B4E9","#E69F00","#D55E00")) +
  labs(title = "Post-Normalization: PCA of Samples",
       subtitle = paste0("PC1 (", var_explained[1], "%) vs PC2 (", var_explained[2], "%)"),
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Save outputs
outdir <- "Results/QC_Post"
ggsave(file.path(outdir, "PCA_Normalized_Data.png"), p_pca, width = 7, height = 5, dpi = 200)
ggsave(file.path(outdir, "PCA_Normalized_Data.pdf"), p_pca, width = 7, height = 5)

# Re-run filtering (in case it didn't save earlier)
X <- exprs(eset_norm)
row_med_abs <- apply(abs(X), 1, median, na.rm = TRUE)
thresh <- quantile(row_med_abs, 0.20, na.rm = TRUE)
min_samples <- ceiling(ncol(X) * 0.20)
keep <- rowSums(abs(X) > thresh) >= min_samples
eset_filt <- eset_norm[keep, ]

# Now count remaining transcripts
transcripts_remaining <- nrow(exprs(eset_filt))
cat("Transcripts remaining after filtering:", transcripts_remaining, "\n")


# Submission Instructions:
# Upload your R script implementing the workflow to GitHub and 
# provide the repository link in the form.

# Google Form: https://forms.gle/1e9tj2Mqf5T9FKEJ7

# Deadline: Sunday 5th October, 2025 (Midnight)

