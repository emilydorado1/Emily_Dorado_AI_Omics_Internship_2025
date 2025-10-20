# =====================================================================
#               AI and Biotechnology / Bioinformatics
#              AI and Omics Research Internship (2025)
#             Module II: Introduction to Genomics Data Analysis
#                        Microarray Data Analysis
# =====================================================================
# From-scratch, turnkey pipeline:
# - Makes fresh folders
# - Downloads GEO data (or optional local upload)
# - Probe -> Gene mapping (auto-annotation)
# - limma differential expression
# - Volcano + Heatmap
# - Saves CSVs + PNGs + results summary
# =====================================================================

# -----------------------
# 0) PROJECT SETUP
# -----------------------
options(stringsAsFactors = FALSE)
set.seed(123)

PROJECT_DIR  <- "Microarray_Assignment_2025"
RESULTS_DIR  <- file.path(PROJECT_DIR, "Results")
PLOTS_DIR    <- file.path(RESULTS_DIR, "Plots")
RDATA_DIR    <- file.path(PROJECT_DIR, "RData")
dir.create(PROJECT_DIR, showWarnings = FALSE)
dir.create(RDATA_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)

# -----------------------
# 1) PARAMETERS
# -----------------------
# Toggle: use a public GEO series (recommended) or upload your own file
USE_LOCAL         <- FALSE     # set TRUE to choose a local file (CSV/TSV)
LOCAL_FILE_PATH   <- NULL      # will be set by file.choose() if USE_LOCAL=TRUE

# A reliable, fully worked example (Gastric cancer vs normal, Affy GPL570)
GSE_ID            <- "GSE79973"

# Labeling for groups (edit if you switch dataset)
CONTROL_LABEL     <- "normal"
CASE_LABEL        <- "cancer"
CONTRAST_NAME     <- paste0(CASE_LABEL, "_vs_", CONTROL_LABEL)

# DEG thresholds
P_CUTOFF          <- 0.05
LOGFC_CUTOFF      <- 1

# Heatmap size
TOP_N_HEATMAP     <- 25

# -----------------------
# 2) PACKAGE INSTALLS
# -----------------------
install_if_missing <- function(pkgs) {
  to_get <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_get)) install.packages(to_get, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc_install_if_missing <- function(pkgs) {
  to_get <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_get)) BiocManager::install(to_get, ask = FALSE, update = FALSE)
}

install_if_missing(c("dplyr","tibble","ggplot2","pheatmap","readr"))
bioc_install_if_missing(c("GEOquery","Biobase","limma","AnnotationDbi"))

library(GEOquery)
library(Biobase)
library(limma)
library(AnnotationDbi)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(readr)

# -----------------------
# 3) DATA INGEST
# -----------------------

exprs_mat <- NULL
pheno_df  <- NULL
platform  <- NULL
annotation_pkg <- NULL

if (!USE_LOCAL) {
  # ---- Pull from GEO (recommended) ----
  message("Downloading GEO series: ", GSE_ID)
  gset_list <- getGEO(GSE_ID, GSEMatrix = TRUE, AnnotGPL = FALSE)
  if (length(gset_list) > 1) {
    # choose the ExpressionSet with the most samples
    ns <- vapply(gset_list, ncol, integer(1))
    eset <- gset_list[[which.max(ns)]]
  } else {
    eset <- gset_list[[1]]
  }
  saveRDS(eset, file.path(RDATA_DIR, paste0(GSE_ID, "_eset.rds")))
  
  # Expression matrix & phenotype
  exprs_mat <- Biobase::exprs(eset)
  pheno_df  <- Biobase::pData(eset)
  platform  <- Biobase::annotation(eset)  # e.g., "GPL570"
} else {
  # ---- Local upload path (expects a *wide* matrix: first col=probe IDs, others=samples) ----
  message("Choose your local expression file (CSV/TSV). First column = Probe IDs; other columns = samples.")
  LOCAL_FILE_PATH <- file.choose()
  # Try to guess delimiter and read
  ext <- tools::file_ext(LOCAL_FILE_PATH)
  if (tolower(ext) %in% c("csv")) {
    raw_df <- readr::read_csv(LOCAL_FILE_PATH, show_col_types = FALSE)
  } else {
    raw_df <- readr::read_tsv(LOCAL_FILE_PATH, show_col_types = FALSE)
  }
  stopifnot(ncol(raw_df) >= 3)     # at least probe + 2 samples
  stopifnot(any(grepl("^probe|^Probe|^PROBE|^ID", colnames(raw_df)[1])))
  
  probe_col <- 1
  rownames(raw_df) <- raw_df[[probe_col]]
  exprs_mat <- as.matrix(raw_df[,-probe_col])
  storage.mode(exprs_mat) <- "double"
  
  # Minimal placeholder phenotype (edit this if using local files)
  pheno_df <- data.frame(
    sample = colnames(exprs_mat),
    group  = rep(c(CONTROL_LABEL, CASE_LABEL), length.out = ncol(exprs_mat)),
    stringsAsFactors = FALSE, row.names = colnames(exprs_mat)
  )
  platform <- NA_character_
}

# -----------------------
# 4) AUTO-INSTALL ANNOTATION PKG (based on platform)
# -----------------------
# Map common GPLs to their AnnotationDbi packages
gpl_to_pkg <- c(
  "GPL570"  = "hgu133plus2.db",  # Affymetrix U133 Plus 2.0
  "GPL96"   = "hgu133a.db",
  "GPL97"   = "hgu133b.db",
  "GPL6244" = "hugene10sttranscriptcluster.db",
  "GPL10558"= "illuminaHumanv4.db",
  "GPL6884" = "illuminaHumanv3.db"
)

if (!is.na(platform) && platform %in% names(gpl_to_pkg)) {
  annotation_pkg <- gpl_to_pkg[[platform]]
  bioc_install_if_missing(annotation_pkg)
  suppressPackageStartupMessages(library(annotation_pkg, character.only = TRUE))
} else if (!is.na(platform)) {
  message("Unknown platform (", platform, "). Trying GPL's annotation table directly...")
} else {
  message("Local file mode: you will be prompted later to provide mapping if needed.")
}

# -----------------------
# 5) PROBE -> GENE MAPPING
# -----------------------
# Strategy:
#  - If we have a known *.db package, use AnnotationDbi::mapIds(PROBEID -> SYMBOL)
#  - Else, try getGEO(GPL) to pull the platform table and map via "ID"->"Gene Symbol"
probe_ids <- rownames(exprs_mat)

map_with_db <- function(pkg, ids) {
  AnnotationDbi::mapIds(
    get(pkg), keys = ids, keytype = "PROBEID",
    column = "SYMBOL", multiVals = "first"
  )
}

gene_symbols <- NULL

if (!is.null(annotation_pkg)) {
  gene_symbols <- map_with_db(annotation_pkg, probe_ids)
} else if (!is.na(platform)) {
  # Fallback: fetch GPL table and map by columns
  gpl <- GEOquery::getGEO(platform, AnnotGPL = TRUE)
  gtab <- Table(gpl)
  # Try common column names for probe and symbol
  probe_col  <- c("ID","ID_REF","ProbeID")
  symbol_col <- c("Gene Symbol","Symbol","SYMBOL","Gene Symbol.1")
  pcol <- probe_col[probe_col %in% colnames(gtab)][1]
  scol <- symbol_col[symbol_col %in% colnames(gtab)][1]
  if (is.na(pcol) || is.na(scol)) {
    stop("Could not find appropriate columns in GPL table for mapping.")
  }
  symbol_map <- gtab[,c(pcol, scol)]
  colnames(symbol_map) <- c("PROBEID","SYMBOL")
  # join
  gene_symbols <- setNames(symbol_map$SYMBOL, symbol_map$PROBEID)[probe_ids]
} else {
  stop("No platform/annotation info available. Provide a mapping file or use a supported GEO series.")
}

# Build annotated data.frame
gene_map_df <- data.frame(
  PROBEID = probe_ids,
  SYMBOL  = unname(gene_symbols),
  stringsAsFactors = FALSE
)

# Remove probes lacking symbols
keep_idx <- !is.na(gene_map_df$SYMBOL) & nchar(gene_map_df$SYMBOL) > 0
exprs_mat <- exprs_mat[keep_idx, , drop = FALSE]
gene_map_df <- gene_map_df[keep_idx, , drop = FALSE]

# Collapse multiple probes per gene by averaging (limma::avereps)
averaged <- limma::avereps(exprs_mat, ID = gene_map_df$SYMBOL)

# -----------------------
# 6) GROUPS / DESIGN
# -----------------------
# Try to auto-detect group labels for GSE79973 (you can edit if using a different series)
if (!USE_LOCAL) {
  # Heuristic: look for a phenotype field with "source_name_ch1" (common in GEO)
  group_vec <- NULL
  if ("source_name_ch1" %in% colnames(pheno_df)) {
    src <- tolower(pheno_df$source_name_ch1)
    if (GSE_ID == "GSE79973") {
      # expected: "gastric mucosa" (normal) vs "gastric adenocarcinoma" (cancer)
      grp <- ifelse(grepl("mucosa", src), CONTROL_LABEL,
                    ifelse(grepl("adenocarcinoma", src), CASE_LABEL, NA))
      group_vec <- grp
    }
  }
  # Fallback to "title" if needed
  if (is.null(group_vec) || any(is.na(group_vec))) {
    tit <- tolower(pheno_df$title)
    grp <- ifelse(grepl("normal|control", tit), CONTROL_LABEL,
                  ifelse(grepl("tumou?r|cancer|case|patient", tit), CASE_LABEL, NA))
    group_vec <- ifelse(is.na(group_vec), grp, group_vec)
  }
  # If still NA, stop with a clear message
  if (any(is.na(group_vec))) {
    stop("Could not auto-detect groups. Edit the group logic for your dataset.")
  }
  names(group_vec) <- rownames(pheno_df)
} else {
  # Local file placeholder groups were created above in pheno_df
  group_vec <- pheno_df$group
  names(group_vec) <- rownames(pheno_df)
}

# Make sure order matches columns of averaged matrix
group_vec <- group_vec[colnames(averaged)]

groups <- factor(group_vec, levels = c(CONTROL_LABEL, CASE_LABEL))
if (nlevels(groups) != 2) stop("Need exactly 2 groups. Check your labels.")

# Design/fit
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
fit <- lmFit(averaged, design)

contrast_matrix <- makeContrasts(contrasts = paste0(CASE_LABEL, "-", CONTROL_LABEL), levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix))

# -----------------------
# 7) DEG TABLES
# -----------------------
deg <- topTable(fit2, number = Inf, adjust.method = "BH")
deg$threshold <- ifelse(deg$adj.P.Val < P_CUTOFF & deg$logFC >  LOGFC_CUTOFF, "Upregulated",
                        ifelse(deg$adj.P.Val < P_CUTOFF & deg$logFC < -LOGFC_CUTOFF, "Downregulated", "No"))

up   <- subset(deg, threshold == "Upregulated")
down <- subset(deg, threshold == "Downregulated")
updown <- rbind(up, down)

# Save CSVs
write.csv(deg,    file.path(RESULTS_DIR, "DEGs_all.csv"), row.names = TRUE)
write.csv(up,     file.path(RESULTS_DIR, "DEGs_up.csv"),  row.names = TRUE)
write.csv(down,   file.path(RESULTS_DIR, "DEGs_down.csv"),row.names = TRUE)
write.csv(updown, file.path(RESULTS_DIR, "DEGs_updown_only.csv"), row.names = TRUE)

# -----------------------
# 8) PLOTS
# -----------------------
# Volcano
volcano_path <- file.path(PLOTS_DIR, "volcano_plot.png")
png(volcano_path, width = 2000, height = 1500, res = 300)
ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No" = "grey70")) +
  theme_minimal() +
  labs(title = paste0("Volcano: ", CONTRAST_NAME),
       x = "log2 Fold Change", y = "-log10(adj. P-value)", color = "Regulation")
dev.off()

# Heatmap (top N by adj.P)
top_genes <- rownames(head(deg[order(deg$adj.P.Val), ], TOP_N_HEATMAP))
hm_mat <- averaged[top_genes, , drop = FALSE]

# Column labels as group_index
grp <- as.character(groups)
colnames(hm_mat) <- ave(grp, grp, FUN = function(x) paste0(x, "_", seq_along(x)))

heatmap_path <- file.path(PLOTS_DIR, paste0("heatmap_top", TOP_N_HEATMAP, ".png"))
png(heatmap_path, width = 2200, height = 1600, res = 300)
pheatmap(hm_mat,
         scale = "row",
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         main = paste0("Top ", TOP_N_HEATMAP, " DEGs (", CONTRAST_NAME, ")"),
         fontsize_row = 6, fontsize_col = 9,
         color = colorRampPalette(c("blue","white","red"))(100))
dev.off()

# -----------------------
# 9) SUMMARY (4–5 lines)
# -----------------------
summary_lines <- c(
  paste0("Dataset: ", ifelse(USE_LOCAL, "Local file", GSE_ID),
         ifelse(!is.na(platform), paste0(" (Platform: ", platform, ")"), "")),
  "Probe→gene mapping: multiple probes per gene were collapsed using limma::avereps (mean of probes sharing the same SYMBOL).",
  paste0("Contrast: ", CONTRAST_NAME, " using limma with empirical Bayes moderation."),
  paste0("Adjusted p<", P_CUTOFF, " & |log2FC|>", LOGFC_CUTOFF, ": ",
         nrow(up), " upregulated, ", nrow(down), " downregulated genes."),
  paste0("Outputs: CSVs in ", RESULTS_DIR, "; plots in ", PLOTS_DIR, ".")
)
writeLines(summary_lines, con = file.path(RESULTS_DIR, "RESULTS_SUMMARY.txt"))

# -----------------------
# 10) README
# -----------------------
readme <- c(
  "# Microarray Assignment (From Scratch)",
  "",
  "Folders:",
  paste0("- Results: ", RESULTS_DIR),
  paste0("- Plots:   ", PLOTS_DIR),
  paste0("- RData:   ", RDATA_DIR),
  "",
  "How to rerun:",
  "1) Open a new RStudio project.",
  "2) Copy this script into the project root and source it.",
  "3) To use a different GEO series, change GSE_ID.",
  "4) To use a local expression matrix, set USE_LOCAL=TRUE; choose file when prompted.",
  "",
  "Notes:",
  "- Expects a two-group comparison. Edit CONTROL_LABEL/CASE_LABEL detection if needed.",
  "- If groups aren’t detected automatically, adjust the logic near section 6."
)
writeLines(readme, file.path(PROJECT_DIR, "README.txt"))

message("Done. Deliverables are in: ", normalizePath(RESULTS_DIR))


cat("Upregulated:", nrow(up), "\nDownregulated:", nrow(down))
